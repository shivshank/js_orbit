'use strict'
// requires THREE js to be defined before loading
// hooray for old school JS and global variables

// all base parameters approximate, based on ISS and Earth (or wherever Google sourced its query
// responses from (since I'm lazy))
// inclination of ISS orbit is not based on anything, I just picked a vector with magnitude 8000

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
        // traveling at ~7,500 m/s
        vel: new THREE.Vector3(0, 0, 1).normalize().multiplyScalar(8500),
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

var KeplerOrbitPrototype = {
    /// given a parameter angle theta find a position on the orbit
    /// theta is an arbitrary parameter (I think!)
    sample_ellipse: function(theta) {
        return this.ellipse.major_vector.clone()
            .multiplyScalar(Math.cos(theta))
            .addScaledVector(this.ellipse.minor_vector, Math.sin(theta))
            .add(this.ellipse.constant);
    },
    apoapsis: function() {
        return this.ellipse.major_vector.clone()
            .multiplyScalar(-1)
            .add(this.ellipse.constant);
    },
    periapsis: function() {
        return this.ellipse.major_vector.clone()
            .add(this.ellipse.constant);
    },
    // this isn't exactly right but oh well. computed radius does not visually line up with anything
    ascending_node: function() {
        // true_anomaly_at_asc
        var v = 2 * Math.PI - this.argument_of_periapsis;
        var r = (this.ang_momentum.lengthSq() / this.mu) / (1 + this.eccentricity * Math.cos(v));
        return this.node_vector.clone().normalize().multiplyScalar(r);
    }
};

/// Find the kepler orbit of b around a.
///
/// assume the orbital plane is the XZ plane
/// for equatorial orbits set the longitude of ascending node to zero
/// all calculations found on wikipedia (look on each parameter's individual page)
///
/// TODO: I believe some of these calculations assume a being at the origin...
/// TODO: Most calculations appear to break down if we play with the parameters. Needs testing.
function calculate_kepler_orbit(a, b) {
    // calculate standard gravitational parameter
    var mu = G * a.mass;
    var pos = b.pos;
    var vel = b.vel;
    var ang_momentum = b.pos.clone().cross(vel);
    // I can't believe this calculation works...:
    // points toward the ascending node from the origin
    var node_vector = (new THREE.Vector3(0, 1, 0)).cross(ang_momentum);

    // this is the vector between periapsis and apoapsis
    var eccentricity_vector = pos.clone().multiplyScalar(vel.lengthSq() / mu - 1 / pos.length())
        .addScaledVector(vel, pos.dot(vel) / mu);

    // our speed is in the thousands of m/s (10^3), so we'll choose our "small number" as 1cm/s
    var eps = 1e-2;
    var longitude_of_ascending_node;
    var argument_of_periapsis;
    if (Math.abs(vel.y) < eps) {
        // equatorial orbit... boo special case!
        longitude_of_ascending_node = 0;
        argument_of_periapsis = Math.atan2(eccentricity_vector.x, eccentricity_vector.z);
        node_vector.x = 0;
        node_vector.y = 0;
        node_vector.z = 1;
    } else {
        // TODO: Is it okay to use ">= zero" for a sign check?
        longitude_of_ascending_node = Math.acos(node_vector.x / node_vector.length());
        if (node_vector.z < 0) {
            longitude_of_ascending_node = 2 * Math.PI - longitude_of_ascending_node;
        }
        argument_of_periapsis = Math.acos(node_vector.dot(eccentricity_vector)
            / node_vector.length() * eccentricity_vector.length());
        if (eccentricity_vector.y < -eps) {
            argument_of_periapsis = 2 * Math.PI - argument_of_periapsis;
        }
    }
    var out = Object.create(KeplerOrbitPrototype);
    var eccentricity = eccentricity_vector.length();
    var inclination = Math.acos(ang_momentum.y / ang_momentum.length());
    var specific_energy = vel.lengthSq() / 2 - mu / pos.length();
    var semi_major_axis = - mu / (2 * specific_energy);

    var longitude_of_periapsis = argument_of_periapsis + longitude_of_ascending_node;
    if (longitude_of_periapsis > 2 * Math.PI) {
        longitude_of_periapsis -= 2 * Math.PI;
    }
    var dist_from_elliptic_center = semi_major_axis * eccentricity;
    // get a vector that points toward periapsis
    // TODO: optimization; my brain isn't working: I think we need to normalize this?
    var major_vector = eccentricity_vector.clone().normalize();
    var constant = major_vector.clone().multiplyScalar(-dist_from_elliptic_center)
        .add(a.pos); 
    major_vector.multiplyScalar(semi_major_axis);
    var minor_vector = major_vector.clone()
        .cross(ang_momentum)
        .normalize()
        .multiplyScalar(semi_major_axis * Math.sqrt(1 - eccentricity * eccentricity));
    Object.assign(out, {
        node_vector,
        eccentricity_vector,
        ang_momentum,
        mu,

        longitude_of_ascending_node,
        argument_of_periapsis,
        eccentricity,
        inclination,
        semi_major_axis,

        ellipse: {
            constant: constant,
            major_vector,
            minor_vector,
        }
    });
    return out;
}