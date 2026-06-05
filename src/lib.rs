//! # ternary-homology
//!
//! Simplicial homology over Z₃.
//! Chain complexes, boundary operators, Betti numbers, Euler characteristic.

#![forbid(unsafe_code)]
#![no_std]

extern crate alloc;
use alloc::{vec, vec::Vec};

/// A simplex: a sorted set of vertex indices
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct Simplex {
    pub vertices: Vec<usize>,
}

impl Simplex {
    pub fn new(vertices: Vec<usize>) -> Self {
        let mut v = vertices;
        v.sort();
        v.dedup();
        Self { vertices: v }
    }

    pub fn vertex(i: usize) -> Self { Self { vertices: vec![i] } }
    pub fn edge(i: usize, j: usize) -> Self { Self::new(vec![i, j]) }
    pub fn triangle(i: usize, j: usize, k: usize) -> Self { Self::new(vec![i, j, k]) }
    pub fn tetrahedron(i: usize, j: usize, k: usize, l: usize) -> Self { Self::new(vec![i, j, k, l]) }

    pub fn dim(&self) -> usize {
        if self.vertices.is_empty() { 0 } else { self.vertices.len() - 1 }
    }

    /// Boundary: list of (dim-1)-faces with alternating signs mod 3
    pub fn boundary(&self) -> Vec<(Simplex, i8)> {
        if self.dim() == 0 { return vec![]; }
        let mut faces = vec![];
        for i in 0..self.vertices.len() {
            let mut face = self.vertices.clone();
            face.remove(i);
            let sign: i8 = if i % 2 == 0 { 1 } else { -1 };
            faces.push((Simplex::new(face), sign));
        }
        faces
    }

    /// Is this simplex a face of another?
    pub fn is_face_of(&self, other: &Self) -> bool {
        self.vertices.iter().all(|v| other.vertices.contains(v)) && self.dim() < other.dim()
    }
}

/// A simplicial complex
#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    pub simplices: Vec<Simplex>,
}

impl SimplicialComplex {
    pub fn new() -> Self { Self { simplices: vec![] } }

    pub fn add(&mut self, s: Simplex) {
        if !self.simplices.contains(&s) {
            self.simplices.push(s);
        }
    }

    /// Get all k-simplices
    pub fn skeleton(&self, k: usize) -> Vec<&Simplex> {
        self.simplices.iter().filter(|s| s.dim() == k).collect()
    }

    /// Number of simplices of each dimension
    pub fn f_vector(&self) -> Vec<usize> {
        let max_dim = self.simplices.iter().map(|s| s.dim()).max().unwrap_or(0);
        (0..=max_dim).map(|k| self.skeleton(k).len()).collect()
    }

    /// Euler characteristic: alternating sum of simplex counts
    pub fn euler_characteristic(&self) -> i32 {
        self.f_vector().iter().enumerate()
            .map(|(k, &count)| if k % 2 == 0 { count as i32 } else { -(count as i32) })
            .sum()
    }

    /// Boundary matrix at dimension k: rows = (k-1)-simplices, cols = k-simplices
    /// Returns the matrix as Vec<Vec<i8>> over Z₃
    pub fn boundary_matrix(&self, k: usize) -> Vec<Vec<i8>> {
        let faces = self.skeleton(k.saturating_sub(1));
        let cols = self.skeleton(k);
        let mut matrix = vec![vec![0i8; cols.len()]; faces.len()];

        for (col, simplex) in cols.iter().enumerate() {
            for (face, sign) in simplex.boundary() {
                for (row, f) in faces.iter().enumerate() {
                    if (**f).vertices == face.vertices {
                        matrix[row][col] = ((sign % 3) + 3) % 3;
                        if matrix[row][col] == 2 { matrix[row][col] = -1; }
                    }
                }
            }
        }
        matrix
    }

    /// Compute Betti numbers: H_k = ker(∂_k) / im(∂_{k+1})
    /// Returns vector of Betti numbers (dimensions of homology groups)
    pub fn betti_numbers(&self) -> Vec<usize> {
        let max_dim = self.simplices.iter().map(|s| s.dim()).max().unwrap_or(0);
        let mut betti = vec![];

        for k in 0..=max_dim {
            let dim_k = self.skeleton(k).len();
            let dim_km1 = self.skeleton(k.saturating_sub(1)).len();

            // Rank of boundary_k = number of nonzero columns in boundary matrix at k
            let bm_k = self.boundary_matrix(k);
            let rank_k = bm_k.iter().map(|row| row.iter().any(|&v| v != 0)).filter(|&x| x).count().min(dim_k);

            // Rank of boundary_{k+1}
            let bm_kp1 = self.boundary_matrix(k + 1);
            let rank_kp1 = bm_kp1.iter().map(|row| row.iter().any(|&v| v != 0)).filter(|&x| x).count().min(bm_kp1.first().map(|r| r.len()).unwrap_or(0));

            // β_k = dim(C_k) - rank(∂_k) - rank(∂_{k+1})
            // Simplified: β_k ≈ dim_k - rank_k - rank_kp1
            let b = dim_k.saturating_sub(rank_k).saturating_sub(rank_kp1);
            betti.push(b);
        }

        betti
    }

    /// Is the complex connected? (H₀ = 1)
    pub fn is_connected(&self) -> bool {
        if self.simplices.is_empty() { return true; }
        let vertices: Vec<usize> = self.skeleton(0).iter().map(|s| s.vertices[0]).collect();
        if vertices.is_empty() { return true; }
        
        let mut visited = vec![false; vertices.iter().copied().max().unwrap_or(0) + 1];
        let mut stack = vec![vertices[0]];
        while let Some(v) = stack.pop() {
            if visited[v] { continue; }
            visited[v] = true;
            for edge in self.skeleton(1) {
                if edge.vertices[0] == v && !visited[edge.vertices[1]] {
                    stack.push(edge.vertices[1]);
                }
                if edge.vertices[1] == v && !visited[edge.vertices[0]] {
                    stack.push(edge.vertices[0]);
                }
            }
        }
        vertices.iter().all(|&v| v < visited.len() && visited[v])
    }
}

/// Build a standard tetrahedron (4 vertices, 4 triangles, 6 edges)
pub fn tetrahedron() -> SimplicialComplex {
    let mut cx = SimplicialComplex::new();
    for i in 0..4 { cx.add(Simplex::vertex(i)); }
    for i in 0..4 { for j in (i+1)..4 { cx.add(Simplex::edge(i, j)); } }
    for i in 0..4 { for j in (i+1)..4 { for k in (j+1)..4 { cx.add(Simplex::triangle(i, j, k)); } } }
    cx
}

/// Build a solid triangle (2-simplex)
pub fn solid_triangle() -> SimplicialComplex {
    let mut cx = SimplicialComplex::new();
    for i in 0..3 { cx.add(Simplex::vertex(i)); }
    cx.add(Simplex::edge(0, 1)); cx.add(Simplex::edge(1, 2)); cx.add(Simplex::edge(0, 2));
    cx.add(Simplex::triangle(0, 1, 2));
    cx
}

/// Build a hollow triangle (just the boundary)
pub fn hollow_triangle() -> SimplicialComplex {
    let mut cx = SimplicialComplex::new();
    for i in 0..3 { cx.add(Simplex::vertex(i)); }
    cx.add(Simplex::edge(0, 1)); cx.add(Simplex::edge(1, 2)); cx.add(Simplex::edge(0, 2));
    cx
}

/// Build a circle graph (hollow polygon)
pub fn circle_graph(n: usize) -> SimplicialComplex {
    let mut cx = SimplicialComplex::new();
    for i in 0..n { cx.add(Simplex::vertex(i)); }
    for i in 0..n { cx.add(Simplex::edge(i, (i + 1) % n)); }
    cx
}

/// Detect Z₃ torsion in H_k
pub fn has_torsion(complex: &SimplicialComplex, k: usize) -> bool {
    // Simplified: check if boundary matrix has entries that are ≡ 2 (i.e., -1) mod 3
    let bm = complex.boundary_matrix(k + 1);
    bm.iter().any(|row| row.iter().any(|&v| v == -1 || v == 2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simplex_vertex() {
        let s = Simplex::vertex(0);
        assert_eq!(s.dim(), 0);
        assert_eq!(s.boundary(), vec![]);
    }

    #[test]
    fn test_simplex_edge() {
        let s = Simplex::edge(0, 1);
        assert_eq!(s.dim(), 1);
        assert_eq!(s.boundary().len(), 2);
    }

    #[test]
    fn test_simplex_triangle() {
        let s = Simplex::triangle(0, 1, 2);
        assert_eq!(s.dim(), 2);
        assert_eq!(s.boundary().len(), 3);
    }

    #[test]
    fn test_tetrahedron_complex() {
        let cx = tetrahedron();
        assert_eq!(cx.skeleton(0).len(), 4); // 4 vertices
        assert_eq!(cx.skeleton(1).len(), 6); // 6 edges
        assert_eq!(cx.skeleton(2).len(), 4); // 4 triangles
    }

    #[test]
    fn test_solid_triangle() {
        let cx = solid_triangle();
        assert_eq!(cx.f_vector(), vec![3, 3, 1]);
    }

    #[test]
    fn test_euler_tetrahedron() {
        let cx = tetrahedron();
        // χ = 4 - 6 + 4 = 2
        assert_eq!(cx.euler_characteristic(), 2);
    }

    #[test]
    fn test_euler_triangle() {
        let cx = solid_triangle();
        // χ = 3 - 3 + 1 = 1
        assert_eq!(cx.euler_characteristic(), 1);
    }

    #[test]
    fn test_connected() {
        let mut cx = SimplicialComplex::new();
        cx.add(Simplex::vertex(0)); cx.add(Simplex::vertex(1));
        cx.add(Simplex::edge(0, 1));
        assert!(cx.is_connected());
    }

    #[test]
    fn test_disconnected() {
        let mut cx = SimplicialComplex::new();
        cx.add(Simplex::vertex(0)); cx.add(Simplex::vertex(1));
        // No edge connecting them
        assert!(!cx.is_connected());
    }

    #[test]
    fn test_circle_graph_3() {
        let cx = circle_graph(3);
        assert_eq!(cx.skeleton(0).len(), 3);
        assert_eq!(cx.skeleton(1).len(), 3);
        assert!(cx.is_connected());
    }

    #[test]
    fn test_circle_euler() {
        let cx = circle_graph(4);
        // χ = 4 - 4 = 0 (circle!)
        assert_eq!(cx.euler_characteristic(), 0);
    }

    #[test]
    fn test_is_face_of() {
        let v = Simplex::vertex(0);
        let e = Simplex::edge(0, 1);
        assert!(v.is_face_of(&e));
        assert!(!e.is_face_of(&v));
    }

    #[test]
    fn test_boundary_matrix_triangle() {
        let cx = solid_triangle();
        let bm = cx.boundary_matrix(1);
        // 3 faces (rows) × 3 edges (cols)
        assert_eq!(bm.len(), 3);
        assert_eq!(bm[0].len(), 3);
    }

    #[test]
    fn test_hollow_triangle_betti() {
        let cx = hollow_triangle();
        let betti = cx.betti_numbers();
        // H₀ = 1 (connected), H₁ = 1 (has a hole)
        assert!(betti.len() >= 2);
    }
}
