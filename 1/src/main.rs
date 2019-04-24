use nalgebra::{DMatrix/*, Dynamic, Matrix, VecStorage*/};
use gnuplot::{Caption, Color, Figure};
use rand::{distributions::Distribution, distributions::Uniform, seq::SliceRandom, Rng};

//type DMatrixi32 = Matrix<i32, Dynamic, Dynamic, VecStorage<i32, Dynamic, Dynamic>>;

fn get_town_rand(size: usize) -> Vec<usize> {
    let mut vec: Vec<usize> = (0..size).collect();
    let mut rng = rand::thread_rng();
    vec.shuffle(&mut rng);
    vec
}

fn get_rand_dist() -> f64 {
    let between = Uniform::new(0, 10000);
    let mut rng = rand::thread_rng();
    let value = between.sample(&mut rng);
    (value / 10000) as f64
}

fn traveling2(d: DMatrix<f64>) -> (f64, Vec<usize>) {
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
        if delta < 0.0 {
            town.swap(c, next1);
            t_dist = t_dist + delta;
            i = 0;
        } else {
            i = i + 1;
        }
    }

    (t_dist, town.to_vec())
}

fn traveling(d: DMatrix<f64>) -> ( f64, Vec<usize>) {

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
    let mut T: f64 = 1.0;

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
        if delta < 0.0 && (-delta as f64 / T).exp() >= get_rand_dist() {
            town.swap(c, next1);
            t_dist = t_dist + delta;
            if delta != 0.0 {
                i = 0;
            }
        } else {
            i = i + 1;
        }

        T = 0.999 * T;
    }

    (t_dist, town.to_vec())
}

/*fn rnd_numb(_i: usize, _j: usize) -> i32 {
    let mut rng = rand::thread_rng();

    rng.gen_range(0, 10)
}*/

fn PARtraveling(n: usize, procs: usize) {
    let x = vec![5.1984, 5.3037, 6.9279, 1.4845, 7.2218];//Vec::new();
    let y = vec![8.5059, 6.6985, 2.4793, 1.5320, 8.5861];//Vec::new();
    let mut array = Vec::new();

    /*for i in 0..n {
        x.push(rnd_numb(0,10));
        y.push(rnd_numb(0,10));
    }*/

    println!("{:?} + {:?}", x, y);

    for i in 0..n {
        for j in 0..n {
            array.push(((x[i] - x[j] as f64).powi(2) + (y[i] - y[j] as f64).powi(2)).sqrt());
        }
    }

    let d = DMatrix::from_row_slice(n, n, &array);

    println!("{:?}", d);

    //let d = DMatrix::from_fn(n, n, rnd_numb);

    let mut xs = Vec::new();
    let mut ys = Vec::new();

    for _ in 0..procs{
        let mut fg = Figure::new();

        let (t_dist, route) = traveling(d.clone());

        println!("{:?}", route);
        println!("{:?}", t_dist);

        for i in 0..n {
            xs.push(x[route[i]]);
            ys.push(y[route[i]]);
        }

        xs.push(x[route[0]]);
        ys.push(y[route[0]]);
        println!("{:?}", xs);
        println!("{:?}", ys);
        
        fg.axes2d()
        .lines(&xs, &ys, &[Caption("Avareage"), Color("black")]);

        fg.show();

        xs.clear();
        ys.clear();
    }

    for _ in 0..procs{
        let mut fg = Figure::new();

        let (_t_dist, route) = traveling2(d.clone());

        for i in 0..n {
            xs.push(x[route[i]]);
            ys.push(y[route[i]]);
        }

        xs.push(x[route[0]]);
        ys.push(y[route[0]]);
        
        fg.axes2d()
        .lines(&xs, &ys, &[Caption("Avareage"), Color("black")]);

        fg.show();

        xs.clear();
        ys.clear();
    }
}

fn main() {
    //let d = DMatrix::from_row_slice(3, 3, &[6, 6, 6, 6, 6, 6, 6, 6, 6]);
    //println!("{:?}", traveling2(d.clone()));
    //println!("{:?}", traveling(d.clone()));

    PARtraveling(5, 2);
}
