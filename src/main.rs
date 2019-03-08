use nalgebra::{DMatrix, Dynamic, Matrix, VecStorage};
use rand::{seq::SliceRandom, Rng};
use std::num;

type DMatrixi32 = Matrix<i32, Dynamic, Dynamic, VecStorage<i32, Dynamic, Dynamic>>;

fn get_town_rand(size: usize) -> Vec<usize> {
    let mut vec: Vec<usize> = (0..size).collect();
    let mut rng = rand::thread_rng();
    vec.shuffle(&mut rng);
    vec
}

fn traveling2(d: DMatrixi32) -> i32 {
    let mut rng = rand::thread_rng();
    let n = d.nrows() - 1;
    let town: &mut [usize] = &mut get_town_rand(n + 1);
    let mut t_dist = d[(town[n], town[0])];

    for i in 0..n {
        t_dist += d[(town[i], town[i + 1])];
    }

    let mut previous: usize;
    let mut next1: usize;
    let mut next2: usize;

    let mut i = 0;
    while i < 100 {
        let c = rng.gen_range(0, n);
        match c {
            _ if c == 0 => {
                previous = n;
                next1 = 1;
                next2 = 2
            }
            _ if c == n - 1 => {
                previous = n - 2;
                next1 = n;
                next2 = 0
            }
            _ if c == n => {
                previous = n - 1;
                next1 = 0;
                next2 = 1
            }
            _ => {
                previous = c - 1;
                next1 = c + 1;
                next2 = c + 2
            }
        };

        let delta = d[(town[previous], town[next1])] + d[(town[c], town[next2])]
            - d[(town[previous], town[c])]
            - d[(town[next1], town[next2])];
        if delta < 0 {
            town.swap(c, next1);
            t_dist = t_dist + delta;
            i = 0;
        } else {
            i = i + 1;
        }
    }

    t_dist
}

fn PARtraveling(n: i32, procs: i32) {
    let x: [i32; n];
    let y: [i32; n];
    let mut rng = rand::thread_rng();

    for i in 0..n {
        x[i] = rng.gen_range(0, 10);
        y[i] = rng.gen_range(0, 10);
    }

    let d: DMatrixi32;

    for i in 0..n {
        for j in 0..n {
            d[(i,j)] = num::sqrt( i32::pow( x[i] - x[j], 2) + i32::pow( y[i] - y[j], 2));
        }
    }

}

fn main() {
    let d = DMatrix::from_row_slice(3, 3, &[6, 6, 6, 6, 6, 6, 6, 6, 6]);
    print!("{:?}", traveling2(d))
}
