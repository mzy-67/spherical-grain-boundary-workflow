import numpy as np
from pymatgen.core.structure import Structure
#from sevenn.calculator import SevenNetCalculator
from pymatgen.core.surface import SlabGenerator
from pymatgen.transformations.site_transformations \
import TranslateSitesTransformation, RemoveSitesTransformation
from pymatgen.transformations.standard_transformations \
import SupercellTransformation
from pymatgen.transformations.standard_transformations import RotationTransformation
from pymatgen.core.lattice import Lattice
from pymatgen.core import Element
from interfacemaster.cellcalc import get_pri_vec_inplane, get_right_hand
from interfacemaster.interface_generator import cross_plane
from interfacemaster.cellcalc import rot
import os
from dataclasses import dataclass, field
from jobflow import Flow, Response, job, Maker
from typing import Callable, Dict, Any, List, Mapping, Optional, Sequence
from pymatgen.io.cif import CifWriter, CifParser
import json
from dataclasses import dataclass
from typing import List, Tuple
import math
import matplotlib.pyplot as plt
from .runtime import run_lammps, spherical_slab_intersection_volume
from .geometry import inclusive_radius_grid, rotation_to_z
#按照最小的位移去移动
def species_symbol(specie) -> str:
    """Return a plain element symbol for a pymatgen specie or string."""
    return getattr(specie, "symbol", str(specie))


def shift_slab_to_origin(slab, anchor_species: Optional[str] = None):
    """Translate the nearest selected atom to the origin.

    ``anchor_species`` is optional.  Leaving it unset makes the geometry code
    material agnostic; setting it (for example to ``Zr`` in the LLZO example)
    reproduces the historical anchoring convention.
    """
    candidates = [site for site in slab if anchor_species is None or species_symbol(site.specie) == anchor_species]
    if not candidates:
        raise ValueError(f"No sites found for origin anchor species {anchor_species!r}")
    coords = np.array([site.coords for site in candidates])
    shift = -coords[np.argmin(np.linalg.norm(coords, axis=1))] + 1e-4
    tt = TranslateSitesTransformation(np.arange(len(slab)), shift, vector_in_frac_coords=False)
    return tt.apply_transformation(slab)
#第一，abc矩阵的读取需要注意，第二输出矩阵 = 原始矩阵的正交归一化版本的转置
def help_matrix(matrix):
    a,b,c = matrix.T
    c_h = np.cross(a,b)
    b_h = np.cross(c_h,a)
    return np.column_stack((a/np.linalg.norm(a), b_h/np.linalg.norm(b_h), c_h/np.linalg.norm(c_h)))
#得倒每一个点的距离
def get_dist_from_mp(structure):
    a, b, c = structure.lattice.matrix
    middle_point = (a+b+c)*1/2
    coords = structure.cart_coords
    return np.linalg.norm(coords - middle_point, axis=1)
#对晶体结构进行对称性分析，并添加一个名为 site_labels的位点属性，该属性包含每个原子的序号#对称等价位置标签和Wyckoff字母
def get_symmetrized_structure(structure):
    return structure.add_site_property('site_labels', np.arange(len(structure)))
#将一个晶体结构绕Z轴旋转指定的角度（以弧度为单位），并返回旋转后的新结构
def get_rotated_structure(structure, R):
    new_lattice = (np.dot(R, structure.lattice.matrix.T)).T
    return Structure(lattice=new_lattice,
                     species=structure.species,
                     coords=structure.frac_coords,
                     site_properties=structure.site_properties,
                    )
#转化坐标系，改成绝对坐标，gpt说将a对应到x，旋转晶格
def align_ab_to_xy(structure):
    old_orient = help_matrix(structure.lattice.matrix.T)
    new_orient= np.eye(3)
    transform = np.linalg.inv(old_orient)
    new_matrix = np.dot(transform, structure.lattice.matrix.T)
    new_lattice = Lattice(new_matrix.T)
    return Structure(lattice=new_lattice,
                     species=structure.species,
                     coords=(np.dot(transform, structure.cart_coords.T)).T,
                     coords_are_cartesian=True,
                    site_properties=structure.site_properties)
#生成一个给定晶体结构和晶面指数（Miller index）的原始表面模型（primitive slab），并将其平移到坐标原点
def get_primitive_slab_hkl(structure, miller_index):
    B  = get_pri_vec_inplane(miller_index,
                             lattice=structure.lattice.matrix.T)
    n = cross_plane(lattice = structure.lattice.matrix.T,
                    n = np.cross(B[:,0], B[:,1]),
                    lim = 10,
                    tol = 0.1,
                    orthogonal = False)
    sup_cart = get_right_hand(np.column_stack((B,n)))
    scaling_matrix = np.dot(np.linalg.inv(structure.lattice.matrix.T), sup_cart).T
    st = SupercellTransformation(scaling_matrix=scaling_matrix)
    return shift_slab_to_origin(align_ab_to_xy(st.apply_transformation(structure)))

def redefine_cubic_lattice(structure):
    return Structure(lattice=np.eye(3) * structure.lattice.a,
                     species=structure.species,
                     coords=structure.cart_coords,
                     site_properties=structure.site_properties,
                     coords_are_cartesian = True
                    )

#将一个晶体结构复制成超晶格，使得超晶格的最小包围盒能够容纳一个半径为R的球体
import numpy as np
def replicate_to_sphere_size(structure, R):
    """
    make a supercell to reach a size enclosing a sphere
    with radius R
    """
    a, b, c = structure.lattice.matrix
    n1 = np.cross(b,c)
    h1 = abs(np.dot(a, n1)/np.linalg.norm(n1))
    n2 = np.cross(a,c)
    h2 = abs(np.dot(b, n2)/np.linalg.norm(n2))
    n3 = np.cross(a,b)
    h3 = abs(np.dot(c, n3)/np.linalg.norm(n3))
    dim1 = int(np.ceil(R*2/h1))
    dim2 = int(np.ceil(R*2/h2))
    dim3 = int(np.ceil(R*2/h3))
    st = SupercellTransformation(((dim1, 0, 0), (0, dim2, 0), (0, 0, dim3)))
    return st.apply_transformation(structure), [dim1, dim2, dim3]
#从一个超晶格结构中提取一个球形区域（半径为r）
def sphere_from_structure(structure, r):
    """
    extract a sphere from a supercell
    """
    distance_from_mp = get_dist_from_mp(structure)
    remove_indices = np.where(distance_from_mp > r)[0]
    rmt = RemoveSitesTransformation(remove_indices)
    structure = rmt.apply_transformation(structure)

    return structure
#固定远距离原子：将距离大于阈值rstar的原子设置为在所有方向上固定
#将两个晶体结构（st1和st2）沿Z轴方向组合成一个新的结构，并在它们之间添加指定的间隙（gap）
def combine_two_structures(st1, st2, gap):
    lattice = st1.lattice
    species = st1.species + st2.species

    middle_point_1 = np.dot([0.5,0.5,0.5], st1.lattice.matrix)
    middle_point_2 = np.dot([0.5,0.5,0.5], st2.lattice.matrix)

    coords = np.append(st1.cart_coords,
                       st2.cart_coords + middle_point_1 - middle_point_2 + [0,0,gap],
                       axis = 0)
    site_properties = {'site_labels': st1.site_properties['site_labels'] + st2.site_properties['site_labels']}
    return Structure(lattice=lattice,
                     species=species,
                     coords=coords,
                     coords_are_cartesian=True,
                    site_properties=site_properties)
#晶体结构（structure）转换到指定边长的立方晶格中，同时保持原子在空间中的实际位置不变（即笛卡尔坐标不变）
def to_cubic_lattice_structure(structure, cubic_a):
    a, b, c = structure.lattice.matrix
    middle_point = (a+b+c)*1/2
    lattice = Lattice.cubic(cubic_a)
    a,b,c = lattice.matrix
    new_middle_point = 1/2*(a+b+c)
    shift = new_middle_point - middle_point

    new_structure = Structure(lattice=lattice,
                     species=structure.species,
                     coords=structure.cart_coords,
                     coords_are_cartesian=True,
                     site_properties=structure.site_properties)
    tt = TranslateSitesTransformation(np.arange(len(new_structure)), new_middle_point-middle_point, False)

    return tt.apply_transformation(new_structure)

def save_site_properties_to_json(unique_sps, indices, filename="site_properties.json"):
    """
    保存为JSON格式，保留完整的数据结构
    """
    # 将numpy数组转换为Python原生类型
    def convert_to_serializable(obj):
        if isinstance(obj, (np.integer, np.int64, np.int32)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, (list, tuple)):
            return [convert_to_serializable(item) for item in obj]
        else:
            return obj
    data = {
        'unique_site_properties': convert_to_serializable(unique_sps),
        'atom_indices_by_type': convert_to_serializable(indices),
        'summary': {
            'total_unique_types': len(unique_sps),
            'total_atoms': sum(len(idx_list) for idx_list in indices)
        }
    }
    with open(filename, 'w') as f:
        json.dump(data, f, indent=2)

def add_surface_site_property(structure, rstar):
    dyn_mtx = np.array([[True, True, True]] * len(structure))
    dists = get_dist_from_mp(structure)
    dyn_mtx[dists > rstar] = [False, False, False]
    return structure.add_site_property('selective_dynamics', dyn_mtx)

class SphereGenerator:
    def __init__(self, structure, anchor_species: Optional[str] = None):
        self.structure = shift_slab_to_origin(
            get_symmetrized_structure(structure), anchor_species=anchor_species
        )

    def get_hemisphere(self, r, dr, dxdydz = [0,0,0], R=np.eye(3), normal = [0,0,1], up = True):
        supercell, dims = replicate_to_sphere_size(self.structure, r + dr)
        supercell = get_rotated_structure(supercell, R)
        dxdydz = [0.5 + dxdydz[0]/dims[0], 0.5 + dxdydz[1]/dims[1], 0.5 + dxdydz[2]/dims[2]]
        tt = TranslateSitesTransformation(np.arange(len(supercell)), dxdydz)
        supercell = tt.apply_transformation(supercell)
        sphere = sphere_from_structure(supercell, r)

        r_vecs = sphere.cart_coords - np.dot(sphere.lattice.matrix.T, [0.5,0.5,0.5])
        if up == True:
            remove_ids = np.where(np.dot(r_vecs, normal) >= 0)[0]
        else:
            remove_ids = np.where(np.dot(r_vecs, normal) < 0)[0]
        rmt = RemoveSitesTransformation(remove_ids)
        return rmt.apply_transformation(sphere)

class SphereGBGenerator:
    def __init__(self, strcture, r, dr, anchor_species: Optional[str] = None, type_map: Optional[Mapping[str, int]] = None):
        self.structure = strcture
        self.r = r
        self.dr = dr
        self.anchor_species = anchor_species
        self.type_map = type_map

    def get_sphere_GB(self, R,
                          normal,
                          dxdydz_1 ,
                          dxdydz_2,
                          gap):
        sphg_1 = SphereGenerator(self.structure, self.anchor_species)
        hemisphere_1 = sphg_1.get_hemisphere(self.r, self.dr, dxdydz_1, np.eye(3), normal, True)
        hemisphere_1 = get_rotated_structure(hemisphere_1, rotation_to_z(normal))
        h_1=to_cubic_lattice_structure(hemisphere_1, 2*(self.r+self.dr))
        lattice_1 = np.array(h_1.lattice.matrix)      # 和你函数里的 lattice 对应
        atoms_1 = np.array(h_1.frac_coords)           # 和你函数里的 atoms 对应
        elements_1 = np.array([species_symbol(sp) for sp in h_1.species])  # 和你函数里的 elements 对应
        write_LAMMPS(lattice=lattice_1, atoms=atoms_1, elements=elements_1, filename="hemisphere_1.data", orthogonal=False, type_map=self.type_map)
        np.savetxt('h1_site_indices', h_1.site_properties['site_labels'], fmt = '%i')
        #hemisphere_1.to_file('1_POSCAR')
        sphg_2 = SphereGenerator(self.structure, self.anchor_species)
        hemisphere_2 = sphg_2.get_hemisphere(self.r, self.dr, dxdydz_2, R, normal,  False)
        hemisphere_2 = get_rotated_structure(hemisphere_2, rotation_to_z(normal))
        #hemisphere_2.to_file('2_POSCAR')
        h_2=to_cubic_lattice_structure(hemisphere_2, 2*(self.r+self.dr))
        lattice_2 = np.array(h_2.lattice.matrix)      # 和你函数里的 lattice 对应
        atoms_2 = np.array(h_2.frac_coords)           # 和你函数里的 atoms 对应
        elements_2 = np.array([species_symbol(sp) for sp in h_2.species])  # 和你函数里的 elements 对应
        write_LAMMPS(lattice=lattice_2, atoms=atoms_2, elements=elements_2, filename="hemisphere_2.data", orthogonal=False, type_map=self.type_map)
        np.savetxt('h2_site_indices', h_2.site_properties['site_labels'], fmt = '%i')
        cb_structure = combine_two_structures(hemisphere_1, hemisphere_2, gap)
        return to_cubic_lattice_structure(cb_structure, 2*(self.r+self.dr))

class SphereGBGenerator_bo:
    def __init__(self, strcture, r, dr, anchor_species: Optional[str] = None):
        self.structure = strcture
        self.r = r
        self.dr = dr
        self.anchor_species = anchor_species

    def get_sphere_GB(self, R,
                          normal,
                          dxdydz_1,
                          dxdydz_2,
                          gap):
        sphg_1 = SphereGenerator(self.structure, self.anchor_species)
        hemisphere_1 = sphg_1.get_hemisphere(self.r, self.dr, dxdydz_1, np.eye(3), normal, True)
        hemisphere_1 = get_rotated_structure(hemisphere_1, rotation_to_z(normal))
        #hemisphere_1.to_file('1_POSCAR')
        sphg_2 = SphereGenerator(self.structure, self.anchor_species)
        hemisphere_2 = sphg_2.get_hemisphere(self.r, self.dr, dxdydz_2, R, normal,  False)
        hemisphere_2 = get_rotated_structure(hemisphere_2, rotation_to_z(normal))
        #hemisphere_2.to_file('2_POSCAR')
        cb_structure = combine_two_structures(hemisphere_1, hemisphere_2, gap)
        return to_cubic_lattice_structure(cb_structure, 2*(self.r+self.dr))

def write_LAMMPS(
        lattice,
        atoms,
        elements,
        filename='lmp_atoms_file',
        orthogonal=False,
        type_map: Optional[Mapping[str, int]] = None):
    """Write a LAMMPS atomic data file with a deterministic species map."""
    elements = np.asarray([species_symbol(element) for element in elements])
    if type_map is None:
        ordered_species = list(dict.fromkeys(elements.tolist()))
        type_map = {symbol: index + 1 for index, symbol in enumerate(ordered_species)}
    else:
        type_map = {str(symbol): int(index) for symbol, index in type_map.items()}
    expected = list(range(1, len(type_map) + 1))
    if sorted(type_map.values()) != expected:
        raise ValueError(f"LAMMPS type IDs must be contiguous and one-based: {type_map}")
    unknown = set(np.unique(elements)) - set(type_map)
    if unknown:
        raise ValueError(f"Species missing from type_map: {sorted(unknown)}")
    # list of elements（按 type 顺序只是为了写文件 header 时好看）
    items = sorted(type_map.items(), key=lambda kv: kv[1])
    element_species = np.array([k for k, v in items])
    element_indices = np.array([v for k, v in items])
    # to Cartesian
    atoms = np.dot(lattice, atoms.T).T
    # 用映射直接生成 species_identifiers
    species_identifiers = np.array([type_map[es] for es in elements]).reshape(1, -1)
    # the atom ID
    IDs = np.arange(len(atoms)).reshape(1, -1) + 1
    # get the final format
    Final_format = np.concatenate((IDs.T, species_identifiers.T), axis=1)
    Final_format = np.concatenate((Final_format, atoms), axis=1)
    # define the box
    xhi, yhi, zhi = lattice[0][0], lattice[1][1], lattice[2][2]
    xlo, ylo, zlo = 0, 0, 0
    xy = lattice[:, 1][0]
    xz = lattice[:, 2][0]
    yz = lattice[:, 2][1]
    with open(filename, 'w', encoding="utf-8") as f:
        f.write(
            '#LAMMPS input file of atoms generated by interface_master. '
            'The elements are: ')
        for ei, es in zip(element_indices, element_species):
            f.write(f'{ei} {es} ')
        f.write(f'\n {len(atoms)} atoms \n \n')
        f.write(f'{len(element_species)} atom types \n \n')
        f.write(f'{xlo:.8f} {xhi:.8f} xlo xhi \n')
        f.write(f'{ylo:.8f} {yhi:.8f} ylo yhi \n')
        f.write(f'{zlo:.8f} {zhi:.8f} zlo zhi \n\n')
        if not orthogonal:
            f.write(f'{xy:.8f} {xz:.8f} {yz:.8f} xy xz yz \n\n')
        f.write('Atoms \n \n')
        np.savetxt(f, Final_format, fmt='%i %i %.16f %.16f %.16f')
import numpy as np

def extract_last_column(traj_file, output_file='atomic_energies.dat'):
    """从LAMMPS轨迹文件中提取最后一列（每原子能量）"""
    energies = []
    atom_ids = []
    with open(traj_file, 'r') as f:
        lines = f.readlines()
    # 查找原子数据开始位置
    for i, line in enumerate(lines):
        if 'ITEM: ATOMS' in line:
            # 找到原子数据开始行
            data_start = i + 1
            # 解析列标题
            columns = line.strip().split()[2:]  # 去掉"ITEM: ATOMS"
            print(f"轨迹文件包含的列: {columns}")
            break
    # 提取每原子能量（最后一列）
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if line and not line.startswith('ITEM:'):  # 跳过空行和标题行
            parts = line.split()
            if len(parts) >= len(columns):  # 确保有足够的数据
                atom_id = int(parts[0])
                energy = float(parts[-1])  # 最后一列
                energies.append(energy)
                atom_ids.append(atom_id)
    energies = np.array(energies)
    return atom_ids, energies

def get_sectional_error_by_r(indices, dists, energies, cut_dists, perfect_energies):
    errors_by_r = []
    sectional_errors_by_r = []
    for i in range(len(indices)):
        errors_by_r.append([])
        for j in cut_dists:
            energies_this_element = energies[indices[i]]
            indices_in_r = np.where(dists[indices[i]] < j)[0]
            errors_by_r[i].append(np.sum((energies_this_element[indices_in_r] - perfect_energies[i])))
        sectional_errors_by_r.append(np.array(errors_by_r[i])/(np.pi*cut_dists**2)*16.02176634)
    total_sectional_errors = np.sum(sectional_errors_by_r, axis=0)
    return errors_by_r, sectional_errors_by_r,total_sectional_errors

def read_lammps_atoms(filename):
    import numpy as np
    with open(filename, "r", encoding="utf-8") as f:
        lines = f.readlines()
    xlo = xhi = ylo = yhi = zlo = zhi = None
    for line in lines:
        if "xlo xhi" in line:
            p = line.split()
            xlo, xhi = float(p[0]), float(p[1])
        elif "ylo yhi" in line:
            p = line.split()
            ylo, yhi = float(p[0]), float(p[1])
        elif "zlo zhi" in line:
            p = line.split()
            zlo, zhi = float(p[0]), float(p[1])
    a = xhi - xlo
    b = yhi - ylo
    c = zhi - zlo
    start_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            start_idx = i + 2
            break
    data = np.loadtxt(lines[start_idx:])  # id type x y z
    ids   = data[:, 0].astype(int)
    types = data[:, 1].astype(int)
    coords = data[:, 2:5]  # x y z，笛卡尔坐标
    return ids, types, coords, (a, b, c)

from jobflow import job


########新增函数功能
from pathlib import Path
SECTION_NAMES = {
    "Masses","Atoms","Velocities","Bonds","Angles","Dihedrals","Impropers",
    "Pair Coeffs","Bond Coeffs","Angle Coeffs","Dihedral Coeffs","Improper Coeffs"
}
def _is_section_title(line: str) -> bool:
    s = line.strip()
    if not s:
        return False
    # 标题行通常是 "Atoms" 或 "Atoms # atomic"
    head = s.split()[0]
    return head in {name.split()[0] for name in SECTION_NAMES}
def _read_section_numeric_lines(lines, start_idx):
    """从 start_idx 开始读，直到空行或下一个 section 标题；返回 (end_idx, numeric_lines)."""
    i = start_idx
    # 跳过紧随标题后的空行（LAMMPS 会忽略标题后的下一行/空行，外部解析时也要兼容） :contentReference[oaicite:1]{index=1}
    while i < len(lines) and lines[i].strip() == "":
        i += 1
    numeric = []
    while i < len(lines):
        s = lines[i].rstrip("\n")
        if s.strip() == "":
            break
        if _is_section_title(s):
            break
        numeric.append(s)
        i += 1
    return i, numeric
def sort_atoms_by_id(infile, outfile):
    lines = Path(infile).read_text(encoding="utf-8", errors="ignore").splitlines(True)
    # 找 Atoms 段标题
    atoms_title = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atoms_title = i
            break
    if atoms_title is None:
        raise ValueError("没找到 Atoms 段")
    # 读 Atoms 数值行
    end_idx, atom_lines = _read_section_numeric_lines(lines, atoms_title + 1)
    # 排序：按第一列 atom id
    def get_id(s):
        s = s.split("#", 1)[0].strip()  # 去掉行尾注释
        return int(s.split()[0])
    atom_lines_sorted = sorted(atom_lines, key=get_id)
    # 写回
    new_lines = lines[:atoms_title+1]  # 含 Atoms 标题行
    # 保留标题后原来的空行（如果你不在意也可不保留）
    i = atoms_title + 1
    while i < len(lines) and lines[i].strip() == "":
        new_lines.append(lines[i]); i += 1
    new_lines += [l if l.endswith("\n") else l + "\n" for l in atom_lines_sorted]
    new_lines += lines[end_idx:]  # Atoms 后面剩余内容原样接上
    Path(outfile).write_text("".join(new_lines), encoding="utf-8")
import numpy as np
def extract_last_column_volume(traj_file, output_file='atomic_energies.dat'):
    """从LAMMPS轨迹文件中提取最后一列（每原子能量）"""
    energies = []
    atom_ids = []
    with open(traj_file, 'r') as f:
        lines = f.readlines()
    # 查找原子数据开始位置
    for i, line in enumerate(lines):
        if 'ITEM: ATOMS' in line:
            # 找到原子数据开始行
            data_start = i + 1
            # 解析列标题
            columns = line.strip().split()[2:]  # 去掉"ITEM: ATOMS"
            print(f"轨迹文件包含的列: {columns}")
            break
    # 提取每原子能量（最后一列）
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if line and not line.startswith('ITEM:'):  # 跳过空行和标题行
            parts = line.split()
            if len(parts) >= len(columns):  # 确保有足够的数据
                atom_id = int(parts[0])
                energy = float(parts[-2])  # 最后一列
                energies.append(energy)
                atom_ids.append(atom_id)
    energies = np.array(energies)
    return atom_ids, energies
def extract_last_column_energy(traj_file, output_file='atomic_energies.dat'):
    """从LAMMPS轨迹文件中提取最后一列（每原子能量）"""
    energies = []
    atom_ids = []
    with open(traj_file, 'r') as f:
        lines = f.readlines()
    # 查找原子数据开始位置
    for i, line in enumerate(lines):
        if 'ITEM: ATOMS' in line:
            # 找到原子数据开始行
            data_start = i + 1
            # 解析列标题
            columns = line.strip().split()[2:]  # 去掉"ITEM: ATOMS"
            print(f"轨迹文件包含的列: {columns}")
            break
    # 提取每原子能量（最后一列）
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if line and not line.startswith('ITEM:'):  # 跳过空行和标题行
            parts = line.split()
            if len(parts) >= len(columns):  # 确保有足够的数据
                atom_id = int(parts[0])
                energy = float(parts[-3])  # 最后一列
                energies.append(energy)
                atom_ids.append(atom_id)
    energies = np.array(energies)
    return atom_ids, energies
from io import StringIO
import numpy as np
def read_lammps_atoms_2(filename, return_image=False):
    with open(filename, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()
    # --- box bounds（正交/三斜都兼容：你的文件 xy xz yz = 0 0 0） ---
    xlo = xhi = ylo = yhi = zlo = zhi = None
    xy = xz = yz = 0.0
    for line in lines:
        if "xlo xhi" in line:
            p = line.split()
            xlo, xhi = float(p[0]), float(p[1])
        elif "ylo yhi" in line:
            p = line.split()
            ylo, yhi = float(p[0]), float(p[1])
        elif "zlo zhi" in line:
            p = line.split()
            zlo, zhi = float(p[0]), float(p[1])
        elif "xy xz yz" in line:
            p = line.split()
            xy, xz, yz = float(p[0]), float(p[1]), float(p[2])
    Lx, Ly, Lz = (xhi - xlo), (yhi - ylo), (zhi - zlo)
    # --- 找到 Atoms 段（别写死 "Atoms # atomic"，更稳） ---
    atoms_header = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            atoms_header = i
            break
    if atoms_header is None:
        raise ValueError("没找到 'Atoms' 段。")
    # Atoms 段数据起始：跳过标题行 + 紧随其后的空行
    i = atoms_header + 1
    while i < len(lines) and lines[i].strip() == "":
        i += 1
    # --- 收集 Atoms 的纯数值行：遇到空行 or 新 section 标题就停止 ---
    atom_lines = []
    while i < len(lines):
        s = lines[i].strip()
        if s == "":
            break
        if s[0].isalpha():          # 例如 Velocities / Bonds / ...
            break
        if s.startswith("#"):
            i += 1
            continue
        # 去掉行尾注释
        s = s.split("#", 1)[0].strip()
        atom_lines.append(s)
        i += 1
    if not atom_lines:
        raise ValueError("Atoms 段为空或解析失败。")
    data = np.loadtxt(StringIO("\n".join(atom_lines)))
    if data.ndim == 1:
        data = data[None, :]
    # 你的文件 Atoms 每行通常是：id type x y z nx ny nz
    ids   = data[:, 0].astype(int)
    types = data[:, 1].astype(int)
    coords = data[:, 2:5].astype(float)
    # 可选：读 image flags（nx ny nz）
    img = None
    if return_image:
        if data.shape[1] >= 8:
            img = data[:, 5:8].astype(int)
        else:
            img = np.zeros((data.shape[0], 3), dtype=int)
    # 返回：盒子长度 + tilt（如果你后续要处理三斜盒子会用到）
    box = (Lx, Ly, Lz, xy, xz, yz)
    return (ids, types, coords, img, box) if return_image else (ids, types, coords, box)
import numpy as np
def get_sectional_error_by_r_energy(indices, dists, energies, cut_dists, perfect_energies):
    errors_by_r = []
    sectional_errors_by_r = []
    for i in range(len(indices)):
        errors_by_r.append([])
        for j in cut_dists:
            energies_this_element = energies[indices[i]]
            indices_in_r = np.where(dists[indices[i]] < j)[0]
            errors_by_r[i].append(np.sum((energies_this_element[indices_in_r] - perfect_energies[i])))
        sectional_errors_by_r.append(np.array(errors_by_r[i])/(np.pi*cut_dists**2)*16.02176634)
    total_sectional_errors = np.sum(sectional_errors_by_r, axis=0)
    return errors_by_r, sectional_errors_by_r,total_sectional_errors
def get_sectional_error_by_r_volume(indices, dists, energies, cut_dists, perfect_energies):
    errors_by_r = []
    sectional_errors_by_r = []
    for i in range(len(indices)):
        errors_by_r.append([])
        for j in cut_dists:
            energies_this_element = energies[indices[i]]
            indices_in_r = np.where(dists[indices[i]] < j)[0]
            errors_by_r[i].append(np.sum((energies_this_element[indices_in_r] - perfect_energies[i])))
        sectional_errors_by_r.append(np.array(errors_by_r[i])/(np.pi*cut_dists**2))
    total_sectional_errors = np.sum(sectional_errors_by_r, axis=0)
    return errors_by_r, sectional_errors_by_r,total_sectional_errors
import re
def get_group_atoms(logfile, group):
    with open(logfile, "r") as f:
        text = f.read()
    m = re.search(rf"(\d+)\s+atoms in group\s+{re.escape(group)}\b", text)
    if not m:
        raise ValueError(f"没找到 group {group}")
    return int(m.group(1))


#define optimizer
from skopt import Optimizer, gp_minimize

def cp_updt_dict(old_dict, up_dict):
    if old_dict == None:
        return up_dict
    new_dict = old_dict.copy()
    new_dict.update(up_dict)
    return new_dict

from skopt.space import Real
from tqdm.auto import tqdm
import shutil

@dataclass
class SphericalGBWorkflowMaker(Maker):
    #BO args
    name: str = 'Sphere GB BO'
    trials: int = 10
    base_estimator: str = 'GP'
    acq_func: str = 'EI'
    acq_optimizer: str = 'lbfgs'
    random_state: int = 42
    metadata: Dict[str, Any] = None
    # Material and potential settings
    crystal_structure: Structure = None
    check_point_file: Optional[str] = None
    bulk_energy_traj_file: str = None
    species_order: Optional[Sequence[str]] = None
    species_masses: Mapping[str, float] = field(default_factory=dict)
    mobile_species: Optional[str] = None
    charge_number: float = 1.0
    origin_anchor_species: Optional[str] = None
    pair_style: str = "deepmd"
    pair_style_args: str = "{checkpoint_file}"
    pair_coeff: str = "* *"
    run_transport: bool = True
    anneal_temperature_K: float = 1200.0
    anneal_equilibration_steps: int = 1000
    anneal_quench_steps: int = 20000
    sphere_R: float = 50
    vaccum_thickness: float = 20
    gb_r: float = 40
    rot_axis: List = None
    rot_angle: float = 0
    normal: List = None
    lammps_executable: str = "lmp"
    mobile_radius_A: float = 39.0
    msd_analysis_radius_A: float = 35.0
    msd_slab_half_thickness_A: float = 15.0
    msd_timestep_ps: float = 0.001
    msd_equilibration_steps: int = 4000
    msd_production_steps: int = 50000
    msd_dump_interval: int = 500
    msd_fit_min_ps: float = 10.0
    msd_fit_max_ps: float = 50.0
    msd_temperatures_K: Tuple[float, ...] = (973.15, 1073.15, 1173.15, 1273.15)

    def resolved_species_order(self) -> Tuple[str, ...]:
        configured = tuple(self.species_order or ())
        if configured:
            order = configured
        else:
            order = tuple(dict.fromkeys(species_symbol(sp) for sp in self.crystal_structure.species))
        present = {species_symbol(sp) for sp in self.crystal_structure.species}
        missing = present.difference(order)
        extra = set(order).difference(present)
        if missing or extra:
            raise ValueError(f"species_order mismatch; missing={sorted(missing)}, extra={sorted(extra)}")
        if self.mobile_species is not None and self.mobile_species not in order:
            raise ValueError(f"mobile_species {self.mobile_species!r} is not in species_order {order}")
        return order

    def resolved_mobile_species(self) -> str:
        if self.mobile_species is None:
            raise ValueError("mobile_species is required for MSD/conductivity analysis")
        return self.mobile_species

    def species_type_map(self) -> Dict[str, int]:
        return {symbol: index + 1 for index, symbol in enumerate(self.resolved_species_order())}

    def lammps_mass_commands(self) -> str:
        lines = []
        for symbol, type_id in self.species_type_map().items():
            mass = float(self.species_masses[symbol]) if symbol in self.species_masses else float(Element(symbol).atomic_mass)
            if not np.isfinite(mass) or mass <= 0:
                raise ValueError(f"Invalid atomic mass for {symbol}: {mass}")
            lines.append(f"mass {type_id} {mass:.10g} # {symbol}")
        return "\n".join(lines)

    def lammps_element_names(self) -> str:
        return " ".join(self.resolved_species_order())

    def lammps_pair_style_command(self) -> str:
        args = self.pair_style_args.format(checkpoint_file=self.check_point_file or "").strip()
        return f"pair_style     {self.pair_style}{(' ' + args) if args else ''}"

    def lammps_pair_coeff_command(self) -> str:
        return f"pair_coeff     {self.pair_coeff}".rstrip()

    def energy_radii(self) -> np.ndarray:
        return inclusive_radius_grid(self.gb_r)

    def lammps_input_static_energy(self):
        lammps_input = f"""
# NANOPARTICLE MELTING
#------------------INITIALIZATION------------------
units           metal
atom_style      atomic
dimension       3
boundary        f f f
read_data       gb.data
{self.lammps_mass_commands()}
#------------------FORCE FIELDS------------------
{self.lammps_pair_style_command()}
{self.lammps_pair_coeff_command()}
neighbor        6.0 bin
neigh_modify   every 10 delay 0 check no
#------------------FIX OUTER ATOMS------------------
# 计算盒子中心
variable        cx equal (xlo+xhi)/2
variable        cy equal (ylo+yhi)/2
variable        cz equal (zlo+zhi)/2
# 定义计算能量的命令
compute energy all pe/atom
compute total_energy all reduce sum c_energy
# 输出设置
thermo_style custom step etotal pe ke temp press vol
comm_modify cutoff 24.00
#------------------ENERGY CALCULATION AND OUTPUT------------------
dump            atom_pe all custom 1 gb_energy.lammpstrj id element type x y z c_energy
dump_modify     atom_pe sort id element {self.lammps_element_names()}
run             0  # 运行0步触发输出
"""
        return lammps_input

    def lammps_input_anneal(self, equilibriate_T, seed, equilibriate_steps, quench_steps):
        lammps_input = f"""
# NANOPARTICLE MELTING
#------------------INITIALIZATION------------------
units           metal
atom_style      atomic
dimension       3
boundary        f f f
read_data       gb_final.data
{self.lammps_mass_commands()}
#------------------FORCE FIELDS------------------
{self.lammps_pair_style_command()}
{self.lammps_pair_coeff_command()}
neighbor        6.0 bin
neigh_modify   every 5 delay 0 check no
#------------------FIX OUTER ATOMS------------------
# 计算盒子中心
variable        cx equal (xlo+xhi)/2
variable        cy equal (ylo+yhi)/2
variable        cz equal (zlo+zhi)/2
variable        zlow equal (zlo+zhi)/2-4
variable        zhigh equal (zlo+zhi)/2+4
variable        zlow1 equal (zlo+zhi)/2-2.5
variable        zhigh1 equal (zlo+zhi)/2+2.5
# 使用region命令创建球形区域来选择原子
region          sphere sphere ${{cx}} ${{cy}} ${{cz}} {self.mobile_radius_A}
region          slab block INF INF INF INF ${{zlow}} ${{zhigh}} units box
region          slab_low block INF INF INF INF INF ${{zlow}} units box
region          slab_up block INF INF INF INF ${{zhigh}} INF units box
region          slab_low1 block INF INF INF INF INF ${{zlow1}} units box
region          slab_up1 block INF INF INF INF ${{zhigh1}} INF units box
# 创建固定原子的组（距离中心>10Å的原子）
group           slab region slab
group           slab_low region slab_low
group           slab_up region slab_up
group           slab_low1 region slab_low1
group           slab_up1 region slab_up1
group           sphere region sphere
group           outer subtract all sphere
group           fixed subtract outer slab
group           unnfixed union slab_up1 slab_low1
group           unfixed subtract all slab_low1 slab_up1

timestep 0.001
# 定义计算能量的命令。20ps
compute energy all pe/atom
compute total_energy unfixed reduce sum c_energy
compute v all voronoi/atom
compute ctemp unfixed temp
velocity all zero linear
fix    freeze unnfixed setforce 0.0 0.0 0.0
min_style      cg
thermo_style custom step etotal pe ke temp press vol
comm_modify cutoff 24.00
thermo      100
thermo_modify flush yes
minimize 0.0 0.01 1000000 10000000
unfix freeze

dump            d_unfixed_pre unnfixed custom 1 unfixed_preanneal.lammpstrj id element type x y z c_energy c_v[1] c_v[2]
dump_modify     d_unfixed_pre sort id element {self.lammps_element_names()}
run             0
undump          d_unfixed_pre

# 固定外层原子
fix             frigid_up slab_up1 rigid single torque * off off off
fix             ffix_low slab_low1 setforce 0.0 0.0 0.0
velocity unfixed create {equilibriate_T} {seed} mom yes dist gaussian
variable q_stp equal "{quench_steps}/1000"
variable e_stp equal "{equilibriate_steps}/1000"
thermo 50
thermo_style custom step c_ctemp c_total_energy
#dump            atom_pe all custom 100 outputs/gb.lammpstrj.* id #element type x y z c_energy
#dump_modify     atom_pe sort id element {self.lammps_element_names()}
#melt
fix f_nvt unfixed nvt temp {equilibriate_T} {equilibriate_T} $(100.0*dt)
run {equilibriate_steps}
unfix f_nvt
#quench
fix f_nvt unfixed nvt temp {equilibriate_T} 1 $(100.0*dt)
run {quench_steps}
unfix f_nvt
unfix frigid_up
unfix ffix_low

dump            d_unfixed_post unnfixed custom 1 unfixed_postanneal.lammpstrj id element type x y z c_energy c_v[1] c_v[2]
dump_modify     d_unfixed_post sort id element {self.lammps_element_names()}
run             0
undump          d_unfixed_post

#zero velocity
velocity all zero linear
fix             freeze fixed setforce 0.0 0.0 0.0
#------------------STRUCTURE OPTIMIZATION------------------
# 设置能量最小化参数
min_style      cg
# 输出设置
thermo_style custom step etotal pe ke temp press vol
comm_modify cutoff 24.00
thermo          100
thermo_modify flush yes
# 在创建新的dump之前，先取消之前的dump
#undump atom_pe
# 执行能量最小化（结构优化）
minimize 0.0 0.01 1000000 10000000
#------------------ENERGY CALCULATION AND OUTPUT------------------
dump            atom_pe all custom 1 gb_energy_anneal.lammpstrj id element type x y z c_energy c_v[1] c_v[2]
dump_modify     atom_pe sort id element {self.lammps_element_names()}
write_data      optimized_gb.data
run             0  # 运行0步触发输出
"""
        return lammps_input

    def lammps_input_surface_energy(self, datafile, trajfile , optdata):
        lammps_input = f"""
# NANOPARTICLE MELTING
#------------------INITIALIZATION------------------
units           metal
atom_style      atomic
dimension       3
boundary        f f f
read_data       {datafile}
{self.lammps_mass_commands()}
#------------------FORCE FIELDS------------------
{self.lammps_pair_style_command()}
{self.lammps_pair_coeff_command()}
neighbor        6.0 bin
neigh_modify   every 5 delay 0 check no
#------------------FIX OUTER ATOMS------------------
# 计算盒子中心
variable        cx equal (xlo+xhi)/2
variable        cy equal (ylo+yhi)/2
variable        cz equal (zlo+zhi)/2
variable        zlow equal (zlo+zhi)/2-4
variable        zhigh equal (zlo+zhi)/2+4
variable        zlow1 equal (zlo+zhi)/2-2.5
variable        zhigh1 equal (zlo+zhi)/2+2.5
# 使用region命令创建球形区域来选择原子
region          sphere sphere ${{cx}} ${{cy}} ${{cz}} {self.mobile_radius_A}
region          slab block INF INF INF INF ${{zlow}} ${{zhigh}} units box
region          slab_low block INF INF INF INF INF ${{zlow}} units box
region          slab_up block INF INF INF INF ${{zhigh}} INF units box
region          slab_low1 block INF INF INF INF INF ${{zlow1}} units box
region          slab_up1 block INF INF INF INF ${{zhigh1}} INF units box
# 创建固定原子的组（距离中心>10Å的原子）
group           slab region slab
group           slab_low region slab_low
group           slab_up region slab_up
group           slab_low1 region slab_low1
group           slab_up1 region slab_up1
group           sphere region sphere
group           outer subtract all sphere
group           fixed subtract outer slab
group           unfixed subtract all slab_low1 slab_up1
timestep 0.001
# 定义计算能量的命令。20ps
compute energy all pe/atom
compute total_energy unfixed reduce sum c_energy
compute v all voronoi/atom
compute ctemp unfixed temp
velocity all zero linear
fix    freeze fixed setforce 0.0 0.0 0.0
min_style      cg
thermo_style custom step etotal pe ke temp press vol
comm_modify cutoff 24.00
thermo      100
thermo_modify flush yes
minimize 0.0 0.01 1000000 10000000
unfix freeze
#------------------ENERGY CALCULATION AND OUTPUT------------------
dump            atom_pe all custom 1 {trajfile} id element type x y z c_energy c_v[1] c_v[2]
dump_modify     atom_pe sort id element {self.lammps_element_names()}
write_data      {optdata}
run             0  # 运行0步触发输出
"""
        return lammps_input

    def lammps_input_MSD(self, T):
        lammps_input = f"""
# -------- 初始化与力场 --------
units           metal
atom_style      atomic
dimension       3
boundary        f f f
read_data       gb_sorted.data
{self.lammps_mass_commands()}
{self.lammps_pair_style_command()}
{self.lammps_pair_coeff_command()}
timestep       {self.msd_timestep_ps}           # ps (metal units)
neighbor       6.0 bin
neigh_modify   every 5 delay 0 check yes
# -------- 定义核心/表面原子并固定表面 --------
# 盒子中心
variable cx equal (xlo+xhi)/2.0
variable cy equal (ylo+yhi)/2.0
variable cz equal (zlo+zhi)/2.0
variable dz equal {self.msd_slab_half_thickness_A}
variable zlow  equal ${{cz}}-${{dz}}
variable zhigh equal ${{cz}}+${{dz}}
variable Rcore equal {self.mobile_radius_A}
variable rcore equal {self.msd_analysis_radius_A}
region  slab block INF INF INF INF ${{zlow}} ${{zhigh}} units box
region  core_outer  sphere ${{cx}} ${{cy}} ${{cz}} ${{Rcore}} units box
region  core_inner  sphere ${{cx}} ${{cy}} ${{cz}} ${{rcore}} units box
group   core   region core_outer          # 可动的“核心原子”
group   core_inner   region core_inner
group   fixed  subtract all core           # 外层壳：固定不动
group   slab  region slab
group   mobile       type {self.species_type_map()[self.resolved_mobile_species()]}
group   mobile_core intersect mobile core_inner slab
group   mobile_all   intersect mobile core
# 把固定层速度清零 + 力清零
velocity fixed set 0.0 0.0 0.0
fix freeze fixed setforce 0.0 0.0 0.0
# Mobile species in the dynamic core
# -------- 热浴参数 & 当前温度 --------
variable Tdamp equal 100.0*dt              # ~100 步的温度驰豫时间
variable T equal {T}
print "===== Run MSD at T = ${{T}} K ====="
# -------- 1. 在该温度下先平衡一段时间 --------
velocity core create ${{T}} 12345 mom yes dist gaussian
fix f_nvt_eq core nvt temp ${{T}} ${{T}} ${{Tdamp}}
thermo       100
thermo_style custom step temp pe
run {self.msd_equilibration_steps}
unfix f_nvt_eq
# -------- 2. 重新计时，打开 MSD 统计 + 生产模拟 --------
reset_timestep 0
dump d_traj all custom {self.msd_dump_interval} traj_T${{T}}.lammpstrj id element type x y z
# 每 {self.msd_dump_interval} 步输出一次所有原子的 id, type, x, y, z
dump_modify d_traj sort id element {self.lammps_element_names()}
# Compute MSD for the configured mobile species in the analysis region
compute msd_mobile mobile_core msd com yes    # c_msd_mobile[1..4]
compute msd_core mobile_all msd com yes
# 每 100 步输出一次瞬时 MSD（不再做时间窗口平均）
fix f_msd all ave/time 1 1 1 c_msd_mobile[1] c_msd_mobile[2] c_msd_mobile[3] c_msd_mobile[4] c_msd_core[1] c_msd_core[2] c_msd_core[3] c_msd_core[4] file msd_T${{T}}.dat
# 输出文件每行： time  MSDx  MSDy  MSDz  MSD_total
# 生产模拟：默认 50000 * 0.001 ps = 50 ps，与论文方法一致
fix f_nvt_prod core nvt temp ${{T}} ${{T}} ${{Tdamp}}
thermo       1000
thermo_style custom step temp c_msd_mobile[4]
run {self.msd_production_steps}
unfix f_nvt_prod
# 关掉 MSD 相关
unfix f_msd
#unfix f_rdf_li_li
#unfix f_rdf_li_o
#unfix f_rdf_li_la
#unfix f_rdf_li_zr
uncompute msd_mobile
uncompute msd_core
#uncompute rdf_li_li
#uncompute rdf_li_o
#uncompute rdf_li_la
#uncompute rdf_li_zr
"""
        return lammps_input

    def sample(self, params):
        #x1, y1, z1, x2, y2, z2, gap = params
        unit,x2, y2, z2, gap = params
        x1,y1,z1=unit*self.norm_unit
        #generate gb
        sgg = SphereGBGenerator_bo(self.crystal_structure, self.sphere_R, self.vaccum_thickness, self.origin_anchor_species)
        gb = sgg.get_sphere_GB(rot(self.rot_axis, self.rot_angle),
                                    self.normal,
                                    [x1,y1,z1],
                                    [x2,y2,z2],
                                    gap)
        ##write lammps structure file
        lattice = np.array(gb.lattice.matrix)      # 和你函数里的 lattice 对应
        atoms = np.array(gb.frac_coords)           # 和你函数里的 atoms 对应
        elements = np.array([species_symbol(sp) for sp in gb.species])  # 和你函数里的 elements 对应
        write_LAMMPS(
            lattice=lattice,
            atoms=atoms,
            elements=elements,
            filename="gb.data",
            orthogonal=False,
            type_map=self.species_type_map(),
        )
        ##gb_indices
        all_indices = np.array(gb.site_properties['site_labels'])
        indices =[]
        for i in range(len(self.crystal_structure)):
            indices.append(np.where(all_indices == i)[0])

        #lammps input file
        with open('lammps.in','w') as f:
            f.write(self.lammps_input_static_energy())

        run_lammps(self.lammps_executable, "lammps.in", "log.lammps")

        #read gb energy
        atom_ids, energies = extract_last_column_energy(self.bulk_energy_traj_file)
        atom_ids1, energies1 = extract_last_column('gb_energy.lammpstrj')

        ids, types, coords, (a, b, c) = read_lammps_atoms("gb.data")

        center = np.array([a/2, b/2, c/2])           # 盒子中心
        dists = np.linalg.norm(coords - center, axis=1)
        errors_by_r, sectional_errors_by_r,total_sectional_errors = \
    get_sectional_error_by_r(indices, dists, energies1, self.energy_radii(), energies)
        shutil.move('gb.data', f'gb_{self.count}.data')
        self.count += 1
        return total_sectional_errors[-1]

    def make(self):
        BO_job = self.BO()
        BO_job.update_metadata({'key':f'{self.metadata}_bys'})
        anneal_job = self.anneal(BO_job.output['best_x7'])
        anneal_job.update_metadata({'key':f'{self.metadata}_anneal'})
        jobs = [BO_job, anneal_job]
        if self.run_transport:
            MSD_job = self.MSD(anneal_job.output['pwd'])
            MSD_job.update_metadata({'key':f'{self.metadata}_MSD'})
            jobs.append(MSD_job)
        return Flow(jobs)

    @job
    def BO(self):
        self.count = 0
        def trial_with_progress(func, n_calls, *args, **kwargs):
            with tqdm(total=n_calls, desc="BO optimizing") as progress:
                def wrapped_func(*call_args, **call_kwargs):
                    result = func(*call_args, **call_kwargs)
                    progress.update(1)
                    return result
                return gp_minimize(
                    wrapped_func,
                    search_space,
                    n_calls=n_calls,
                    n_random_starts=max(1, int(0.1 * n_calls)),
                    base_estimator=self.base_estimator,
                    acq_func=self.acq_func,
                    acq_optimizer=self.acq_optimizer,
                    *args,
                    **kwargs,
                )
        #norm = np.linalg.norm(self.normal)
        self.norm_unit = self.normal / np.linalg.norm(self.normal)
        search_space = [Real(0, 0.5*self.crystal_structure.lattice.a, name = 'unit'),
                 Real(0, 0.5*self.crystal_structure.lattice.a, name = 'x2'),
                 Real(0, 0.5*self.crystal_structure.lattice.b, name = 'y2'),
                 Real(0, 0.5*self.crystal_structure.lattice.c, name = 'z2'),
                 Real(0, 3, name = 'gap')]

        result = trial_with_progress(self.sample, n_calls=self.trials, random_state=self.random_state)
        x7_iters = []
        for unit, x2, y2, z2, gap in result.x_iters:
            x1, y1, z1 = (unit * self.norm_unit).tolist()
            x7_iters.append([x1, y1, z1, x2, y2, z2, gap])
        best_unit, best_x2, best_y2, best_z2, best_gap = result.x
        best_x1, best_y1, best_z1 = (best_unit * self.norm_unit).tolist()
        best_x7 = [best_x1, best_y1, best_z1, best_x2, best_y2, best_z2, best_gap]

        return {
             "x": x7_iters,          # 你要的 7 参数轨迹
             "y": result.func_vals.tolist(),
             "best_x7": best_x7 ,
         "pwd": os.getcwd()}
       # return {'x':result.x_iters, 'y':result.func_vals}

    @job
    def anneal(self, cs):
        sgg = SphereGBGenerator(self.crystal_structure, self.sphere_R, self.vaccum_thickness, self.origin_anchor_species, self.species_type_map())
        gb = sgg.get_sphere_GB(rot(self.rot_axis, self.rot_angle),
                                    self.normal,
                                    [cs[0],cs[1],cs[2]],
                                    [cs[3],cs[4],cs[5]],
                                    cs[6])
        lattice = np.array(gb.lattice.matrix)      # 和你函数里的 lattice 对应
        atoms = np.array(gb.frac_coords)           # 和你函数里的 atoms 对应
        elements = np.array([species_symbol(sp) for sp in gb.species])  # 和你函数里的 elements 对应
        write_LAMMPS(
            lattice=lattice,
            atoms=atoms,
            elements=elements,
            filename="gb_final.data",
            orthogonal=False,
            type_map=self.species_type_map(),
        )
        np.savetxt('gb_site_indices', gb.site_properties['site_labels'], fmt = '%i')
        ##gb_indices
        #lammps input file
        with open('lammps_anneal.in','w') as f:
            f.write(self.lammps_input_anneal(equilibriate_T=self.anneal_temperature_K, seed=12345, equilibriate_steps=self.anneal_equilibration_steps, quench_steps=self.anneal_quench_steps))
        run_lammps(self.lammps_executable, "lammps_anneal.in", "log_anneal.lammps")

        for h_id in [1,2]:
            with open(f'lammps_h{h_id}.in','w') as f:
                f.write(self.lammps_input_surface_energy(datafile = f'hemisphere_{h_id}.data',
                                                        trajfile = f'h{h_id}_energy.lammpstrj',
                                                        optdata=f'optimized_h{h_id}.data'))
            run_lammps(
                self.lammps_executable,
                f"lammps_h{h_id}.in",
                f"log_h{h_id}.lammps",
            )
        sort_atoms_by_id("optimized_gb.data", "gb_sorted.data")
        sort_atoms_by_id("optimized_h1.data", "h1_sorted.data")
        sort_atoms_by_id("optimized_h2.data", "h2_sorted.data")
        indices = []
        for i in range(len(self.crystal_structure)):
            indice = []
            with open("gb_site_indices", "r", encoding="utf-8") as f:
                for line_no, line in enumerate(f, start=0):  # 行号从 1 开始
                    line = line.strip()
                    if not line:
                        continue           # 跳过空行
                    x = int(line)          # 如果是小数就用 float(line)
                    if x == i:
                        indice.append(line_no)
            indices.append(indice)
        indices1 = []
        for i in range(len(self.crystal_structure)):
            indice1 = []
            with open("h1_site_indices", "r", encoding="utf-8") as f:
                for line_no, line in enumerate(f, start=0):  # 行号从 1 开始
                    line = line.strip()
                    if not line:
                        continue           # 跳过空行
                    x = int(line)          # 如果是小数就用 float(line)
                    if x == i:
                        indice1.append(line_no)
            indices1.append(indice1)
        indices2 = []
        for i in range(len(self.crystal_structure)):
            indice2 = []
            with open("h2_site_indices", "r", encoding="utf-8") as f:
                for line_no, line in enumerate(f, start=0):  # 行号从 1 开始
                    line = line.strip()
                    if not line:
                        continue           # 跳过空行
                    x = int(line)          # 如果是小数就用 float(line)
                    if x == i:
                        indice2.append(line_no)
            indices2.append(indice2)

        ids, types, coords, box = read_lammps_atoms_2("gb_sorted.data")
        center = np.array([box[0]/2, box[1]/2, box[2]/2])           # 盒子中心
        dists = np.linalg.norm(coords - center, axis=1)
        ids, types, coords, box = read_lammps_atoms_2("h1_sorted.data")
        center = np.array([box[0]/2, box[1]/2, box[2]/2])           # 盒子中心
        dists1 = np.linalg.norm(coords - center, axis=1)
        ids, types, coords, box = read_lammps_atoms_2("h2_sorted.data")
        center = np.array([box[0]/2, box[1]/2, box[2]/2])           # 盒子中心
        dists2 = np.linalg.norm(coords - center, axis=1)
        #read gb energy
        atom_ids, energies = extract_last_column_energy(r'gb_energy_anneal.lammpstrj')
        atom_ids1, energies1 = extract_last_column_energy(r'h1_energy.lammpstrj')
        atom_ids2, energies2 = extract_last_column_energy(r'h2_energy.lammpstrj')
        atom_ids0, energies0= extract_last_column_energy(self.bulk_energy_traj_file)
        atom_ids, volume0 = extract_last_column_volume(self.bulk_energy_traj_file)
        atom_ids1, volume = extract_last_column_volume(r'gb_energy_anneal.lammpstrj')
        errors_by_r, sectional_errors_by_r,total_sectional_errors = get_sectional_error_by_r_energy(indices, dists, energies, self.energy_radii(), energies0)
        errors_by_r1, sectional_errors_by_r1,total_sectional_errors1 = get_sectional_error_by_r_energy(indices1, dists1, energies1, self.energy_radii(), energies0)
        errors_by_r2, sectional_errors_by_r2,total_sectional_errors2 = get_sectional_error_by_r_energy(indices2, dists2, energies2, self.energy_radii(), energies0)
        total_sectional_errors_all=np.array(total_sectional_errors1)+np.array(total_sectional_errors2)-np.array(total_sectional_errors)
        errors_by_r_volume, sectional_errors_by_r_vomume,total_sectional_errors_volume = get_sectional_error_by_r_volume(indices, dists, volume, self.energy_radii(), volume0)
        selected_index = len(self.energy_radii()) - 1
        return {
             "radii_A": self.energy_radii().tolist(),
             "gb_energy_profile": np.asarray(total_sectional_errors).tolist(),
             "work_of_separation_profile": np.asarray(total_sectional_errors_all).tolist(),
             "excess_volume_profile": np.asarray(total_sectional_errors_volume).tolist(),
             "selected_radius_A": float(self.energy_radii()[selected_index]),
             "gb_energy": float(total_sectional_errors[selected_index]),
             "work_of_separation": float(total_sectional_errors_all[selected_index]),
             "excess_volume": float(total_sectional_errors_volume[selected_index]),
             # Legacy aliases for downstream notebooks.
             "wb": np.asarray(total_sectional_errors_all).tolist(),
             "excess": np.asarray(total_sectional_errors_volume).tolist(),
             "pwd": os.getcwd() }
    @job
    def MSD(self,dir_):
        import shutil
        dir = Path(dir_)
        src = f"{dir_}/gb_sorted.data"   # 原文件完整路径
        dst = "./gb_sorted.data"                  # 目标路径（当前目录）
        shutil.copy(src, dst)
        import numpy as np
        for T in self.msd_temperatures_K:
            with open(f'lammps_{T}.in','w') as f:
                f.write(self.lammps_input_MSD(T))
            run_lammps(
                self.lammps_executable,
                f"lammps_{T}.in",
                f"log_{T}.lammps",
            )

        x = get_group_atoms(f"log_{self.msd_temperatures_K[0]}.lammps", "mobile_core")

        MSD_FILES = [
            {"filename": f"msd_T{T}.dat", "T_K": T}
            for T in self.msd_temperatures_K
        ]
        # Number density uses the intersection of the configured sphere and centered slab.
        analysis_volume_A3 = spherical_slab_intersection_volume(
            self.msd_analysis_radius_A,
            self.msd_slab_half_thickness_A,
        )
        results, extra = run_msd_analysis(
            MSD_FILES,
            n_mobile=x,
            volume_A3=analysis_volume_A3,
            time_factor_ps=self.msd_timestep_ps,
            fit_t_min_ps=self.msd_fit_min_ps,
            fit_t_max_ps=self.msd_fit_max_ps,
            charge_number=self.charge_number,
        )
        T_25C = 298.15
        sigma_25_S_m = None
        a = None
        b = None
        if extra["arrhenius_slope"] is not None:
            sigma_25_S_m = sigma_from_arrhenius(
                T_25C, extra["n_m3"], extra["arrhenius_slope"], extra["arrhenius_intercept"], self.charge_number
            )
            sigma_25_S_cm = sigma_25_S_m * 0.01
            #print(f"预测 25°C (T={T_25C:.2f} K) 电导率：")
            #print(f"  σ(298.15K) = {sigma_25_S_m:.3e} S/m  ≈ {sigma_25_S_cm:.3e} S/cm")

            # 同时把公式打印出来（可选）
            a = extra["arrhenius_slope"]
            b = extra["arrhenius_intercept"]
            #print("\nσ(T) 公式（基于拟合 ln D = a*(1/T)+b）：")
            #print("  D(T) = exp(b + a/T)")
            #print("  σ(T) = (n*e^2/(kB*T)) * exp(b + a/T)")
            #print(f"  其中 a={a:.6e}, b={b:.6e}")
        else:
            print("\n温度点不足（至少2个）无法外推 25°C 电导率。")

        return {
             "x": [MSDResult_to_dict(dct) for dct in results],
            "y": extra,
          "pwd": os.getcwd(),
          "Tc":sigma_25_S_m,
          "a":a,
         "b":b}

@dataclass
class MSDResult:
    T_K: float
    D_SI: float          # m^2/s  (用于严格物理计算)
    D_cm2_s: float       # cm^2/s (方便输出 & 画图)
    D_raw: float         # 原始单位：Å^2/ps
    sigma_S_m: float     # S/m
    sigma_S_cm: float    # S/cm


K_B_J_PER_K = 1.380649e-23
E_CHARGE_C = 1.602176634e-19


# ================== 核心函数 ==================
def load_msd_file(
    filename: str,
    time_col: int = 0,
    msd_col: int = 4,
    time_factor_ps: float = 0.001,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    从 MSD 文件里读取时间（ps）和 MSD（Å^2）
    假设：注释行以 '#' 开头
    """
    data = np.loadtxt(filename, comments='#')
    if len(data.shape) == 1:
        time_raw = np.array([data[time_col]])
        msd_raw = np.array([data[msd_col]])
    else:
        time_raw = data[:, time_col]
        msd_raw = data[:, msd_col]

    time_ps = time_raw * time_factor_ps
    msd_A2  = msd_raw                     # 这里假设本来就是 Å^2
    return time_ps, msd_A2
def fit_D_from_MSD(time_ps: np.ndarray, msd_A2: np.ndarray,
                t_min_ps: float = None, t_max_ps: float = None,
                dim: int = 3) -> float:
    """
    用 Einstein 关系从 MSD 拟合 D（返回单位：Å^2/ps）
    MSD(t) ~ 2*d*D*t  => 斜率 m = 2*d*D  => D = m / (2*d)
    """
    # 选择拟合区间
    if t_min_ps is None and t_max_ps is None:
        # 默认后半段数据
        n = len(time_ps)
        idx_fit = np.arange(n // 2, n)
    else:
        idx_fit = np.where(
            (time_ps >= (t_min_ps if t_min_ps is not None else time_ps[0])) &
            (time_ps <= (t_max_ps if t_max_ps is not None else time_ps[-1]))
        )[0]

    t_fit = time_ps[idx_fit]
    msd_fit = msd_A2[idx_fit]

    if len(t_fit) < 2:
        raise ValueError(
            f"MSD fit interval contains {len(t_fit)} point(s); at least 2 are required"
        )

    # 简单线性拟合：MSD = m * t + b
    #print(t_fit, msd_fit)
    coeffs = np.polyfit(t_fit, msd_fit, 1)
    slope = coeffs[0]     # Å^2 / ps
    # 爱因斯坦关系
    D_raw = slope / (2.0 * dim)  # Å^2/ps

    return D_raw
def D_raw_to_SI(D_raw_A2_per_ps: float) -> float:
    """
    把 D 从 Å^2/ps 转成 m^2/s
    1 Å = 1e-10 m; 1 Å^2 = 1e-20 m^2; 1 ps = 1e-12 s
    => 1 Å^2/ps = 1e-8 m^2/s
    """
    return D_raw_A2_per_ps * 1e-8
def D_raw_to_cm2s(D_raw_A2_per_ps: float) -> float:
    """
    把 D 从 Å^2/ps 转成 cm^2/s
    1 Å = 1e-8 cm; 1 Å^2 = 1e-16 cm^2; 1 ps = 1e-12 s
    => 1 Å^2/ps = 1e-4 cm^2/s
    """
    return D_raw_A2_per_ps * 1e-4
def compute_number_density(n_mobile: int, volume_A3: float) -> float:
    """
    计算载流子数密度 n (1/m^3)
    内部仍然用 SI 单位；输出时可以同时给出 cm^-3 方便你看
    """
    if n_mobile <= 0:
        raise ValueError("n_mobile must be positive")
    if volume_A3 <= 0:
        raise ValueError("volume_A3 must be positive")
    V_cell_m3 = volume_A3 * 1e-30
    n = n_mobile / V_cell_m3  # 1/m^3
    return n
def nernst_einstein_sigma(D_SI: float, n: float, T_K: float, charge_number: float = 1.0) -> float:
    """
    Nernst-Einstein 关系：
    σ = n q^2 D / (k_B T)
    返回单位：S/m
    """
    return n * ((abs(charge_number) * E_CHARGE_C) ** 2) * D_SI / (K_B_J_PER_K * T_K)
def fit_Ea_from_D(results: List[MSDResult]) -> Tuple[float, float]:
    """
    用 Arrhenius 关系从 D(T) 拟合迁移能垒 Ea：
    D = D0 * exp(-Ea / (k_B T))
    => ln D = ln D0 - Ea/(k_B) * 1/T
    拟合 ln D vs 1/T 的斜率，Ea = -slope * k_B
    返回：Ea_J (J), Ea_eV (eV)
    """
    T_list = np.array([r.T_K for r in results])
    D_list = np.array([r.D_SI for r in results])  # 用 SI 或 cm^2/s 拟合都可以，只差一个常数
    x = 1.0 / T_list             # 1/K
    y = np.log(D_list)           # ln(D)
    slope, intercept = np.polyfit(x, y, 1)
    Ea_J = -slope * K_B_J_PER_K
    Ea_eV = Ea_J / E_CHARGE_C
    return Ea_J, Ea_eV, slope, intercept
def D_arrhenius_SI(T_K: float, slope: float, intercept: float) -> float:
    """
    由 Arrhenius 拟合参数恢复 D(T)（单位：m^2/s）
    ln D = slope*(1/T) + intercept
    """
    return float(np.exp(intercept + slope / T_K))

def sigma_from_arrhenius(T_K: float, n: float, slope: float, intercept: float, charge_number: float = 1.0) -> float:
    """
    用 Arrhenius 的 D(T) + Nernst-Einstein 得到 σ(T)（单位：S/m）
    """
    D_SI = D_arrhenius_SI(T_K, slope, intercept)
    return nernst_einstein_sigma(D_SI, n, T_K, charge_number=charge_number)
# ================== 主流程 ==================
def run_msd_analysis(
    MSD_FILES,
    *,
    n_mobile: int,
    volume_A3: float,
    time_factor_ps: float,
    fit_t_min_ps: float,
    fit_t_max_ps: float,
    dimension: int = 3,
    time_col: int = 0,
    msd_col: int = 4,
    charge_number: float = 1.0,
) -> Tuple[List[MSDResult], Dict[str, Any]]:
    """
    返回：
    - results: 每个温度一个 MSDResult
    - extra: 额外信息（如 n、Ea 拟合参数等）
    """
    n = compute_number_density(n_mobile, volume_A3)
    n_cm3 = n / 1e6
    results = []
    for item in MSD_FILES:
        filename = item["filename"]
        T_K = item["T_K"]
        time_ps, msd_A2 = load_msd_file(
            filename,
            time_col=time_col,
            msd_col=msd_col,
            time_factor_ps=time_factor_ps,
        )
        #print(time_ps, msd_A2)
        D_raw = fit_D_from_MSD(
            time_ps, msd_A2,
            t_min_ps=fit_t_min_ps,
            t_max_ps=fit_t_max_ps,
            dim=dimension,
        )

        if not np.isfinite(D_raw) or D_raw <= 0:
            raise ValueError(f"Non-positive diffusion coefficient fitted from {filename}: {D_raw}")

        D_SI = D_raw_to_SI(D_raw)
        D_cm2_s = D_raw_to_cm2s(D_raw)
        sigma_S_m = nernst_einstein_sigma(D_SI, n, T_K, charge_number=charge_number)
        sigma_S_cm = sigma_S_m * 0.01
        results.append(MSDResult(
            T_K = T_K,
            D_SI = D_SI,
            D_cm2_s = D_cm2_s,
             D_raw = D_raw,
            sigma_S_m = sigma_S_m,
            sigma_S_cm = sigma_S_cm)
            )

    extra: Dict[str, Any] = {
        "n_m3": n,
        "n_cm3": n_cm3,
        "n_mobile": n_mobile,
        "analysis_volume_A3": volume_A3,
        "charge_number": charge_number,
        "fit_t_min_ps": fit_t_min_ps,
        "fit_t_max_ps": fit_t_max_ps,
        "Ea_J": None,
        "Ea_eV": None,
        "arrhenius_slope": None,
        "arrhenius_intercept": None,
    }

    if len(results) >= 2:
        Ea_J, Ea_eV, slope, intercept = fit_Ea_from_D(results)
        extra.update({
            "Ea_J": Ea_J,
            "Ea_eV": Ea_eV,
            "arrhenius_slope": slope,
            "arrhenius_intercept": intercept,
        })
    return results, extra

def MSDResult_to_dict(msd_result):
    return {'T_K':msd_result.T_K,
            'D_SI':msd_result.D_SI,
           'D_cm2_s':msd_result.D_cm2_s,
           'D_raw':msd_result.D_raw,
           'sigma_S_m':msd_result.sigma_S_m,
           'sigma_S_cm':msd_result.sigma_S_cm}


# Backward-compatible name used by the original LLZO release.
SpheregbBOMaker = SphericalGBWorkflowMaker
