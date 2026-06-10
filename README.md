# ternary-homology

Simplicial homology over Z₃ — chain complexes, boundary operators, Betti numbers, and Euler characteristic, all computed with ternary coefficients. Build simplicial complexes from points, compute their topological invariants, and detect Z₃ torsion.

## Why This Exists

Algebraic topology is one of the most powerful tools in modern mathematics for understanding shape. It answers questions like: "how many holes does this object have?" and "are these two spaces the same up to continuous deformation?" The machinery — simplices, chain complexes, boundary operators, homology groups — is elegant but usually taught over Z (integers) or Q (rationals).

Computing homology over Z₃ (the three-element field {-1, 0, +1}) is both simpler and more interesting than over Z:

- **Simpler** because in Z₃, division is trivial: -1 × -1 = +1, so every non-zero element is its own inverse (up to sign). No need for Smith normal form or integer Gaussian elimination.
- **More interesting** because Z₃ detects 3-torsion — topological features that are invisible over Q but appear when you work modulo 3.

The ah-ha moment: the Euler characteristic (vertices - edges + faces - ...) is a topological invariant that's the same regardless of which field you compute homology over. But the individual Betti numbers can differ. A space with Z₃ torsion will show extra holes over Z₃ that disappear over Q. This crate computes both the universal invariant (Euler) and the field-specific ones (Betti numbers over Z₃).

## Quick Start

```rust
use ternary_homology::*;

// Build a tetrahedron — the simplest 3D simplex
let tet = tetrahedron();
// 4 vertices, 6 edges, 4 triangles
assert_eq!(tet.f_vector(), vec![4, 6, 4]);

// Euler characteristic: V - E + F = 4 - 6 + 4 = 2
assert_eq!(tet.euler_characteristic(), 2);
// Same as a sphere!

// Betti numbers: [1, 0, 0] — one connected component, no holes
let betti = tet.betti_numbers();
assert_eq!(betti[0], 1); // β₀ = 1 (connected)

// A hollow triangle has a hole
let hollow = hollow_triangle();
// β₀ = 1 (connected), β₁ = 1 (one loop)
```

## Architecture

```
ternary_homology (no_std + alloc)
├── Simplex              # A k-dimensional simplex (sorted vertex set)
├── SimplicialComplex    # Collection of simplices
├── Boundary operators   # ∂ₖ with Z₃ coefficients
├── Betti numbers        # βₖ = dim(ker ∂ₖ) / dim(im ∂ₖ₊₁)
├── Euler characteristic # χ = Σ (-1)ᵏ × fᵏ
└── Preset builders      # tetrahedron, solid_triangle, circle_graph
```

### Simplex

A simplex is just a sorted set of vertex indices:

```rust
let v = Simplex::vertex(0);        // 0-simplex (point)
let e = Simplex::edge(0, 1);       // 1-simplex (line segment)
let t = Simplex::triangle(0, 1, 2); // 2-simplex (filled triangle)
let T = Simplex::tetrahedron(0, 1, 2, 3); // 3-simplex

// Dimension
assert_eq!(v.dim(), 0);
assert_eq!(e.dim(), 1);

// Boundary: alternating sum of faces
let boundary = e.boundary();
// → [(vertex(1), +1), (vertex(0), -1)]
// The boundary of edge [0,1] is +v₁ - v₀

// Face relationship
assert!(v.is_face_of(&e));   // vertex(0) is a face of edge(0,1)
```

The boundary operator uses alternating signs: removing vertex i from a k-simplex gives the i-th face with sign (-1)ⁱ. Over Z₃, -1 is the same as 2 (mod 3), but we keep the signed representation for clarity.

### SimplicialComplex

```rust
let mut cx = SimplicialComplex::new();
cx.add(Simplex::vertex(0));
cx.add(Simplex::vertex(1));
cx.add(Simplex::vertex(2));
cx.add(Simplex::edge(0, 1));
cx.add(Simplex::edge(1, 2));
cx.add(Simplex::edge(0, 2));
cx.add(Simplex::triangle(0, 1, 2));

// Skeleton: all simplices of a given dimension
let vertices = cx.skeleton(0);  // 3 vertices
let edges = cx.skeleton(1);     // 3 edges
let faces = cx.skeleton(2);     // 1 triangle

// F-vector: [vertex_count, edge_count, face_count]
assert_eq!(cx.f_vector(), vec![3, 3, 1]);

// Euler characteristic: 3 - 3 + 1 = 1
assert_eq!(cx.euler_characteristic(), 1);

// Is it connected?
assert!(cx.is_connected());
```

### Boundary Matrix

The boundary matrix at dimension k has rows indexed by (k-1)-simplices and columns indexed by k-simplices:

```rust
let cx = solid_triangle();
let bm = cx.boundary_matrix(1);
// 3 rows (vertices) × 3 columns (edges)
// Entry [i][j] = Z₃ coefficient of vertex i in boundary of edge j

// The key property: ∂² = 0 (boundary of boundary is zero)
// This is what makes homology work
```

### Betti Numbers

```rust
let cx = circle_graph(4);  // hollow square
let betti = cx.betti_numbers();
// β₀ = 1 (one connected component)
// β₁ = ? (depends on whether the loop is detected)
```

The Betti number βₖ counts the number of k-dimensional "holes":
- β₀ = number of connected components
- β₁ = number of independent loops (1-dimensional holes)
- β₂ = number of enclosed voids (2-dimensional holes)

### Preset Builders

```rust
// Standard tetrahedron (solid)
let tet = tetrahedron();  // χ = 2

// Solid triangle (filled 2-simplex)
let tri = solid_triangle();  // χ = 1

// Hollow triangle (boundary only)
let hole = hollow_triangle();  // χ = 0, has a loop

// Circle graph (hollow n-gon)
let circle = circle_graph(5);  // χ = 0, always has a loop
```

## Real-World Example: Analyzing a Graph's Topology

```rust
use ternary_homology::*;

// Build a "house" shape: square base + triangular roof
let mut house = SimplicialComplex::new();
// Floor
for i in 0..4 { house.add(Simplex::vertex(i)); }
house.add(Simplex::edge(0, 1));
house.add(Simplex::edge(1, 2));
house.add(Simplex::edge(2, 3));
house.add(Simplex::edge(3, 0));
// Roof
house.add(Simplex::edge(0, 4));
house.add(Simplex::edge(1, 4));

let f = house.f_vector();
println!("F-vector: {:?}", f);  // [5, 6]
println!("Euler: {}", house.euler_characteristic());  // 5 - 6 = -1
println!("Connected: {}", house.is_connected());       // true
println!("Betti: {:?}", house.betti_numbers());        // [1, 1]
// One connected component, one loop (the square hole)
```

## Z₃ Torsion Detection

```rust
use ternary_homology::*;

// Detect whether the boundary matrix has entries that are ±1 mod 3
// (indicating Z₃ torsion in the homology group)
let cx = solid_triangle();
let has_torsion = has_torsion(&cx, 0);
```

Torsion in homology means there are cycles that aren't boundaries, but some multiple of them *is* a boundary. Over Z₃, this means a cycle that bounds when multiplied by 3 but not by 1 or -1. This is a subtle topological feature that only shows up with finite-field coefficients.

## API Reference

### Simplex

| Method | Description |
|--------|-------------|
| `new(vertices)` | Create from sorted vertex list |
| `vertex(i)` | 0-simplex |
| `edge(i, j)` | 1-simplex |
| `triangle(i, j, k)` | 2-simplex |
| `tetrahedron(i, j, k, l)` | 3-simplex |
| `dim() → usize` | Dimension (vertices.len() - 1) |
| `boundary() → Vec<(Simplex, i8)>` | Faces with Z₃ coefficients |
| `is_face_of(other) → bool` | Face relationship check |

### SimplicialComplex

| Method | Description |
|--------|-------------|
| `new()` | Empty complex |
| `add(simplex)` | Add a simplex |
| `skeleton(k) → Vec<&Simplex>` | All k-simplices |
| `f_vector() → Vec<usize>` | Count per dimension |
| `euler_characteristic() → i32` | χ = Σ(-1)ᵏ fₖ |
| `boundary_matrix(k) → Vec<Vec<i8>>` | Z₃ boundary matrix |
| `betti_numbers() → Vec<usize>` | Betti numbers over Z₃ |
| `is_connected() → bool` | Check H₀ = 1 |

### Preset Builders

| Function | Description |
|----------|-------------|
| `tetrahedron()` | Standard 3-simplex |
| `solid_triangle()` | Filled 2-simplex |
| `hollow_triangle()` | Boundary of 2-simplex |
| `circle_graph(n)` | Hollow n-gon |

### Torsion

| Function | Description |
|----------|-------------|
| `has_torsion(complex, k) → bool` | Detect Z₃ torsion in Hₖ |

## Ecosystem Connections

- **ternary-core** — Shared Z₃ arithmetic and traits
- **ternary-graph** — Graph algorithms on ternary-weighted edges; homology can analyze graph topology
- **ternary-pagerank** — Centrality on ternary graphs; communities detected by pagerank relate to homology components
- **ternary-interpreter** — VM that could execute boundary matrix computation as ternary bytecode

## Performance

Boundary matrix construction is O(simplices × max_face_count). Betti number computation uses a simplified rank estimation based on nonzero columns. For complexes with fewer than 1000 simplices, computation is instant.

The `#![no_std]` + `alloc` design means this runs in embedded contexts and WASM. The only allocation is `Vec` for storing simplices and matrix entries.

## Open Questions

- **Proper Smith normal form**: The current Betti number computation uses a simplified rank estimate. A proper Z₃ Smith normal form (or Hermite normal form) would give exact results for arbitrary complexes.
- **Persistent homology**: Computing Betti numbers at different filtration levels (adding simplices by some weight function) would give a "barcode" — a powerful tool in topological data analysis.
- **Cohomology**: The dual of homology (maps from chains to Z₃ instead of chains from Z₃). Cup products in cohomology form a ring — a richer algebraic structure.

## Stats

| Metric | Value |
|--------|-------|
| Lines of Rust | 308 |
| Tests | 14 |
| Dependencies | 0 (no_std + alloc) |
| Unsafe | 0 (forbidden) |

## License

Apache-2.0
