pub struct UnionFind {
    parent: Vec<usize>,
    rank: Vec<usize>,
    size: Vec<usize>,
}

impl UnionFind {
    pub fn new(n: usize) -> UnionFind {
        UnionFind {
            parent: (0..n).map(|i| i).collect::<Vec<usize>>(),
            rank: vec![0; n],
            size: vec![1; n],
        }
    }

    pub fn root(&mut self, x: usize) -> usize {
        if self.parent[x] == x {
            x
        } else {
            self.parent[x] = self.root(self.parent[x]);
            self.parent[x]
        }
    }
    
    pub fn unite(&mut self, x: usize, y: usize) {
        let x = self.root(x);
        let y = self.root(y);
        if x == y {
            return;
        } else {
            if self.rank[x] < self.rank[y] {
                self.parent[x] = y;
                self.size[y] += self.size[x];
            } else {
                self.parent[y] = x;
                self.size[x] += self.size[y];
                if self.rank[x] == self.rank[y] {
                    self.rank[x] += 1;
                }
            }
        }
    }
    
    pub fn same(&mut self, x: usize, y: usize) -> bool {
        self.root(x) == self.root(y)
    }
    
    pub fn size(&mut self, x: usize) -> usize {
        let x_root = self.root(x);
        self.size[x_root]
    }
}


#[derive(Debug, Copy, Clone)]
pub struct Edge {
    pub u: usize,
    pub v: usize,
    pub cost: f64,
    pub sample_num: i64,
}

impl Edge {
    pub fn new(u: usize, v: usize, cost: f64) -> Edge {
        Edge {
            u,
            v,
            cost,
            sample_num: 0,
        }
    }
}

pub fn kruskal(n: usize, es: &mut Vec<Edge>) -> f64 {
    es.sort_by(|s, t| s.cost.partial_cmp(&t.cost).unwrap()); 
    let mut uf = UnionFind::new(n);
    let mut res: f64  = 0.0;
    for e in es {
        if !uf.same(e.u, e.v) {
            uf.unite(e.u, e.v);
            res += e.cost;
        }
    }
    res
}

// return the optimal base I using kruskal's algorithm
pub fn optimal_base(n: usize, es: &mut Vec<(usize, Edge)>) -> Vec<usize> {
    es.sort_by(|(_, s), (_, t)| s.cost.partial_cmp(&t.cost).unwrap()); 
    let mut uf = UnionFind::new(n);
    let mut ret = Vec::new();
    for (i, e) in es {
        if !uf.same(e.u, e.v) {
            uf.unite(e.u, e.v);
            ret.push(i.clone());
        }
    }
    ret
}


// use std::io::Read;
// // AOJ GRL_2_A
// fn main() {
//     let mut buffer = String::new();
//     std::io::stdin().read_to_string(&mut buffer).unwrap();
//     let mut iter = buffer.split_whitespace();
// 
//     let V: usize = iter.next().unwrap().parse().unwrap();
//     let E: usize = iter.next().unwrap().parse().unwrap();
// 
//     let mut edge: Vec<Edge> = (0..E).map(|_| Edge::new(
//         iter.next().unwrap().parse().unwrap(),
//         iter.next().unwrap().parse().unwrap(),
//         iter.next().unwrap().parse().unwrap())
//     ).collect();
//     
//     let res = kruskal(V, &mut edge);
//     println!{"{}", res};
// }