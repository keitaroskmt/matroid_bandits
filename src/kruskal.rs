use std::collections::btree_set::Union;

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


#[derive(Debug)]
pub struct Edge {
    u: usize,
    v: usize,
    cost: i64,
}

impl Edge {
    pub fn new(u: usize, v: usize, cost: i64) -> Edge {
        Edge {
            u,
            v,
            cost,
        }
    }
}

pub fn kruskal(n: usize, es: &mut Vec<Edge>) -> i64 {
    es.sort_by(|s, t| t.cost.cmp(&s.cost));
    let mut uf = UnionFind::new(n);
    let mut res: i64  = 0;
    for e in es {
        if !uf.same(e.u, e.v) {
            uf.unite(e.u, e.v);
            res += e.cost;
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn es_sort() {
        let mut es = Vec::new();
        es.push(Edge::new(0, 1, 2));
        es.push(Edge::new(1, 2, 4));
        es.push(Edge::new(0, 2, 3));
        println!("{:?}", es);
        es.sort_by(|s, t| t.cost.cmp(&s.cost)); 
        println!("{:?}", es);
    }
}
