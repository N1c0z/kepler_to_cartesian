use std::f64::consts::PI;
use std::fs;
use std::f64;
use std::fs::File;
use std::io::ErrorKind;
use std::io::Write;

const AU: f64 = 149597870700f64;
const MU: f64 = 1.32712440041e20;
const STEP_SIZE: usize = 30; //days
const AMOUNT_OF_DAYS: usize = 49590;
const STEP_AMOUNT: usize = AMOUNT_OF_DAYS/STEP_SIZE;

#[derive(Clone, Copy)]
struct Cartesian {
    x: f64,
    y: f64,
    z: f64,
}

#[derive(Clone, Copy)]
struct KeplerElements {
    time: f64,
    index: u8,
    semi_major_axis: f64,
    eccentricity: f64,
    argument_periaps: f64,
    longitude_ascending_node: f64,
    eccentric_anomaly:f64,
    inclination: f64,
    true_anomaly: f64,
    mean_anomaly: f64,
    epoch_t0: f64,
    epoch_t: f64,
}

fn main() {
    let folder_dir = "/Users/nico/Documents/Code stuff/kepler_to_cartesian/aeitest_000000000000.dat".to_string();
    let mut folder_dir2 = "/Users/nico/Documents/Code stuff/kepler_to_cartesian/".to_string();
    let mut celestial_bodies = deserialize(&folder_dir);
    
    let mut position:Vec<[Cartesian;STEP_AMOUNT]> = Vec::new();

    for x in 0..celestial_bodies.len() {
        position.push([Cartesian{x: 0f64, y: 0f64, z: 0f64};STEP_AMOUNT]);
        for y in 0..STEP_AMOUNT {
            let a = calc_algo(celestial_bodies.get_mut(x).unwrap());
            position.get_mut(x).unwrap().get_mut(y).unwrap().x = a.x;
            position.get_mut(x).unwrap().get_mut(y).unwrap().y = a.y;
            position.get_mut(x).unwrap().get_mut(y).unwrap().z = a.z;
            if celestial_bodies.get(x).unwrap().epoch_t == 0f64 {
                celestial_bodies.get_mut(x).unwrap().epoch_t = STEP_SIZE as f64;
            }
        }

        folder_dir2.push_str(&x.to_string());
        folder_dir2.push_str(".txt");
        
        let mut f = File::create(&folder_dir2).unwrap();

        for z in 0..STEP_AMOUNT {
            f.write_all(format!("{}, {}, {}, {},\n",(z * STEP_SIZE) as f64/365.25f64 + 2000f64, position.get(x).unwrap().get(z).unwrap().x/AU, position.get(x).unwrap().get(z).unwrap().y/AU, position.get(x).unwrap().get(z).unwrap().z/AU).as_bytes()).unwrap();
        }
        
        for _ in 0..5 {
            folder_dir2.pop();
        }
    }

}


fn calc_algo(body: &mut KeplerElements) -> Cartesian {
    let delta_t = 86400f64 * (body.epoch_t - body.epoch_t0);
    body.mean_anomaly = body.mean_anomaly + delta_t * (MU/ body.semi_major_axis.powi(3)).sqrt();
    if body.mean_anomaly >= 2f64*PI {
        body.mean_anomaly -= 2f64*PI;
    } 

    body.eccentric_anomaly = newton_raphson(body.mean_anomaly, body.eccentricity, body.mean_anomaly, 50);

    body.true_anomaly = 2f64 * arctan2((1f64+body.eccentricity).sqrt() * (body.eccentric_anomaly/2f64).sin(), (1f64-body.eccentricity).sqrt() * (body.eccentric_anomaly/2f64).cos());

    let distance_to_central_body: f64 = body.semi_major_axis*(1f64 - body.eccentricity * body.eccentric_anomaly.cos());

    let ot: Cartesian = Cartesian { x: distance_to_central_body * body.true_anomaly.cos(), y: distance_to_central_body * body.true_anomaly.sin(), z: 0f64 };

    let mut rt = Cartesian {x: 0f64, y: 0f64, z: 0f64};

    rt.x = ot.x * (body.argument_periaps.cos() * body.longitude_ascending_node.cos() - body.argument_periaps.sin()*body.inclination.cos()*body.longitude_ascending_node.sin()) - ot.y * (body.argument_periaps.sin()*body.longitude_ascending_node.cos() + body.argument_periaps.cos()* body.inclination.cos() * body.longitude_ascending_node.sin());
    rt.y = ot.x * (body.argument_periaps.cos() * body.longitude_ascending_node.sin() + body.argument_periaps.sin()*body.inclination.cos()*body.longitude_ascending_node.cos()) + ot.y * (body.argument_periaps.cos()*body.inclination.cos()*body.longitude_ascending_node.cos()-body.argument_periaps.sin()*body.longitude_ascending_node.sin());
    rt.z = ot.x * body.argument_periaps.sin()*body.inclination.sin()+ot.y * body.argument_periaps.cos()*body.inclination.sin();
    rt
}

fn newton_raphson(eccentric_anomaly0: f64, eccentricity: f64, mean_anomaly: f64, step_amount: u8) -> f64{
    let mut f: f64 = 0f64;
    let mut d_f: f64 = 0f64;
    let mut eccentric_anomaly = eccentric_anomaly0;

    for x in 0..step_amount {
        f = eccentric_anomaly - mean_anomaly - (eccentricity * eccentric_anomaly.sin());
        d_f = 1f64 - eccentricity * eccentric_anomaly.cos();
        eccentric_anomaly = eccentric_anomaly - (f/d_f);
    }
    eccentric_anomaly
}

pub fn arctan2(y: f64, x: f64) -> f64 {
    if x > 0f64 {
        return (y/x).atan();
    }
    if y >= 0f64 && x < 0f64 {
        return (y/x).atan() + PI;
    }
    if y < 0f64 && x <0f64 {
        return (y/x).atan() - PI;
    }
    if y > 0f64 && x == 0f64 {
        return PI/2f64;
    }
    if y < 0f64 && x == 0f64 {
        return PI/2f64;
    }
    panic!("x and y were both 0");
}

fn deserialize<'a>(folder_dir: &'a String) -> Vec<KeplerElements>{

    let content = fs::read_to_string(folder_dir).unwrap();
    let lines: Vec<&str> = content.split('\n').collect();

    //let mut foo = [KeplerElements {eccentric_anomaly: 0f64, true_anomaly: 0f64, time: 0f64, index: 0, semi_major_axis: 0f64, eccentricity: 0f64, argument_periaps: 0f64, longitude_ascending_node: 0f64, inclination: 0f64, mean_anomaly:0f64, epoch_t: 0f64, epoch_t0: 0f64}; 6];
    let mut foo: Vec<KeplerElements> = Vec::new();

    for x in 0..(lines.len()-1) {
        let mut word = lines.get(x).unwrap().split_ascii_whitespace();
        foo.push(KeplerElements {eccentric_anomaly: 0f64, true_anomaly: 0f64, time: 0f64, index: 0, semi_major_axis: 0f64, eccentricity: 0f64, argument_periaps: 0f64, longitude_ascending_node: 0f64, inclination: 0f64, mean_anomaly:0f64, epoch_t: 0f64, epoch_t0: 0f64});
        foo.get_mut(x).unwrap().time = word.next().unwrap().parse::<f64>().unwrap();
        foo.get_mut(x).unwrap().index = word.next().unwrap().parse::<u8>().unwrap();
        //yes
        foo.get_mut(x).unwrap().semi_major_axis = word.next().unwrap().parse::<f64>().unwrap() * 149597870700f64;
        //yes
        foo.get_mut(x).unwrap().eccentricity = word.next().unwrap().parse::<f64>().unwrap();
        //yes
        foo.get_mut(x).unwrap().inclination = word.next().unwrap().parse::<f64>().unwrap();
        //yes
        foo.get_mut(x).unwrap().longitude_ascending_node = word.next().unwrap().parse::<f64>().unwrap();
        //yes
        foo.get_mut(x).unwrap().argument_periaps = word.next().unwrap().parse::<f64>().unwrap();
        
        //foo.get_mut(x).unwrap().true_anomaly = word.next().unwrap().parse::<f64>().unwrap();
        word.next();

        //foo.get_mut(x).unwrap().eccentric_anomaly = word.next().unwrap().parse::<f64>().unwrap();
        word.next();

        //yes
        foo.get_mut(x).unwrap().mean_anomaly = word.next().unwrap().parse::<f64>().unwrap();

    }
    drop(folder_dir);
    foo
}