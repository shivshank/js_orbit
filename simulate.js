'use strict'
// requires THREE js to be defined before loading
// hooray for old school JS and global variables

// all base parameters approximate, based on ISS and Earth (or wherever Google sourced its query
// responses from (since I'm lazy))

var model = {
    planet: {
        // this is the mass of the earth in kg
        mass: 5.972e24,
        pos: new THREE.Vector3(0, 0, 0),
        radius: 6000000
    },
    sat: {
        // mass of ISS
        mass: 400000,
        // ISS is roughly ~400km away from a planet that has radius 6000km
        pos: new THREE.Vector3(6410000, 0, 0),
        // traveling at 7,500 m/s
        vel: new THREE.Vector3(0, 0, 8000),
    }
};

var G = 6.67408e-11;
/// calculate the force of gravity on b between a and b
function f_g(a, b) {
    var delta = a.pos.clone();
    delta.sub(b.pos);
    var r = delta.length();
    delta.normalize()
        .multiplyScalar(G * a.mass * b.mass / (r * r));
    return delta;
}

function step_model(dt) {
    // time_scale is how many seconds to simulate per real life second
    var time_scale = 120;
    var remaining = dt * time_scale;
    dt = time_scale * 1 / 120;
    while (remaining > dt) {
        // we'll just assume the planet is always more massive than the sat
        // such that the force on the planet can be neglected

        // velocity verlet:
        // x(t + dt) = x(t) + v(t) dt + 0.5 a(t) dt^2
        // calculate a(t + dt), assume it does not depend on v(t + dt)
        // v(t + dt) = v(t) + 0.5 (a(t) + a(t + dt)) dt
        var initial_acc_term = f_g(model.planet, model.sat)
            .multiplyScalar(1 / model.sat.mass)
            .multiplyScalar(0.5 * dt * dt);
        var pos = model.sat.pos;
        var vel = model.sat.vel;
        // pos += vel * dt
        pos.add(vel.clone().multiplyScalar(dt));
        // pos += initial_acc
        pos.add(initial_acc_term); 
        // now we have that  pos = x(t + dt) and vel = v(t)
        var second_acc_term = f_g(model.planet, model.sat)
            .multiplyScalar(1 / model.sat.mass)
            .multiplyScalar(0.5 * dt);
        // vel += (0.5 * acc * dt * dt)/dt i.e. vel += 0.5 * acc * dt
        // probably introduces more error than necessary but **** these OOP vectors
        // in a language that discourages allocation and prevents operator overloading >:(
        vel.add(initial_acc_term.multiplyScalar(1 / dt));
        // vel += a(t + dt)
        vel.add(second_acc_term);
        // consume the simulated time
        remaining -= dt;
    }
    return remaining / time_scale;
}