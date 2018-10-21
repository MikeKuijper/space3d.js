class vector {
  constructor(x, y, z) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.projected = 0;

    this.projection = [
      [1, 0, 0],
      [0, 1, 0]
    ];
  }

  project() {
    this.projected = matrixToVector(matrixMultiplyVector(this.projection, this));
    this.projected.z = 0;
  }

  rotate(xAngle, yAngle, zAngle) {
    if (xAngle != 0) this.rotateX(xAngle);
    if (yAngle != 0) this.rotateRelativeY(yAngle);
    if (zAngle != 0) this.rotateRelativeZ(zAngle);
  }

  rotateAround(vector, xAngle, yAngle, zAngle) {
    if (xAngle != 0) this.rotateX(xAngle);
    if (yAngle != 0) this.rotateRelativeY(yAngle);
    if (zAngle != 0) this.rotateRelativeZ(zAngle);
  }

  rotateRelativeX(angle) {
    let rotationX = [
      [1, 0, 0],
      [0, cos(angle), -sin(angle)],
      [0, sin(angle), cos(angle)]
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationX, this.projected));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateRelativeY(angle) {
    let rotationY = [
      [cos(angle), 0, sin(angle)],
      [0, 1, 0],
      [sin(angle), 0, cos(angle)],
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationY, this.projected));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateRelativeZ(angle) {
    let rotationZ = [
      [cos(angle), -sin(angle), 0],
      [sin(angle), cos(angle), 0],
      [0, 0, 1]
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationZ, this.projected));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateX(angle) {
    let rotationX = [
      [1, 0, 0],
      [0, cos(angle), -sin(angle)],
      [0, sin(angle), cos(angle)]
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationX, this));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateY(angle) {
    let rotationY = [
      [cos(angle), 0, sin(angle)],
      [0, 1, 0],
      [sin(angle), 0, cos(angle)],
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationY, this));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateZ(angle) {
    let rotationZ = [
      [cos(angle), -sin(angle), 0],
      [sin(angle), cos(angle), 0],
      [0, 0, 1]
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationZ, this));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateAround(vector, xAngle, yAngle, zAngle) {
    if (xAngle != 0) this.rotateXAround(vector, xAngle);
    if (yAngle != 0) this.rotateRelativeYAround(vector, yAngle);
    if (zAngle != 0) this.rotateRelativeZAround(vector, zAngle);
  }

  rotateRelativeXAround(vector, angle) {
    let rotationX = [
      [1, 0, 0],
      [0, cos(angle), -sin(angle)],
      [0, sin(angle), cos(angle)]
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationX, new vector(this.projected.x - vector.x, this.projected.y - vector.y, this.projected.y - vector.z)));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateRelativeYAround(vector, angle) {
    let rotationY = [
      [cos(angle), 0, sin(angle)],
      [0, 1, 0],
      [sin(angle), 0, cos(angle)],
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationY, new vector(this.projected.x - vector.x, this.projected.y - vector.y, this.projected.y - vector.z)));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateRelativeZAround(vector, angle) {
    let rotationZ = [
      [cos(angle), -sin(angle), 0],
      [sin(angle), cos(angle), 0],
      [0, 0, 1]
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationZ, new vector(this.projected.x - vector.x, this.projected.y - vector.y, this.projected.y - vector.z)));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateXAround(vector, angle) {
    let rotationX = [
      [1, 0, 0],
      [0, cos(angle), -sin(angle)],
      [0, sin(angle), cos(angle)]
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationX, new vector(this.x - vector.x, this.y - vector.y, this.y - vector.z)));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateYAround(vector, angle) {
    let rotationY = [
      [cos(angle), 0, sin(angle)],
      [0, 1, 0],
      [sin(angle), 0, cos(angle)],
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationY, new vector(this.x - vector.x, this.y - vector.y, this.y - vector.z)));
    this.projected = new vector(v.x, v.y, v.z);
  }

  rotateZAround(vector, angle) {
    let rotationZ = [
      [cos(angle), -sin(angle), 0],
      [sin(angle), cos(angle), 0],
      [0, 0, 1]
    ];
    let v = matrixToVector(matrixMultiplyVector(rotationZ, new vector(this.x - vector.x, this.y - vector.y, this.y - vector.z)));
    this.projected = new vector(v.x, v.y, v.z);
  }
}

function vectorToMatrix(v) {
  let m = [
    [v.x],
    [v.y],
    [v.z]
  ];
  return m;
}

function matrixToVector(m) {
  let v = new vector();
  v.x = m[0][0];
  v.y = m[1][0];
  if (m.length > 2) {
    v.z = m[2][0];
  }
  return v;
}

function matrixMultiplyMatrix(a, b) {
  let columnsA = a[0].length;
  let rowsA = a.length;
  let columnsB = b[0].length;
  let rowsB = b.length;

  if (columnsA != rowsB) {
    throw("Number of columns of A must be equal to the number of rows of B");
  }

  let result = createArray(rowsA, columnsB);

  for (let i = 0; i < rowsA; i++) {
    for (let j = 0; j < columnsB; j++) {
      let sum = 0;
      for (let k = 0; k < columnsA; k++) {
        sum += a[i][k] * b[k][j];
      }
      result[i][j] = sum;
    }
  }

  return result;
}

function matrixMultiplyVector(a, b) {
  let m = vectorToMatrix(b);
  return matrixMultiplyMatrix(a, m);
}

function point3d(x, y, z) {
  let v = new vector(x, y, z);
  v.project();
  point(v.projected.x / zoom, v.projected.y / zoom);
}

function ellipse3d(x, y, z, xRotation, yRotation, zRotation, d1, d2) {
  let v = new vector(x / zoom, y / zoom, z / zoom);
  v.rotate(xRotation, yRotation, zRotation);
  ellipse(v.projected.x, v.projected.y, d1 / zoom, d2 / zoom);
}

function text3d(message, x, y, z) {
  let v = new vector(x / zoom, y / zoom, z / zoom);
  v.project();
  text(message, v.projected.x, v.projected.y);
}

function createArray(length) {
  var arr = new Array(length || 0),
    i = length;

  if (arguments.length > 1) {
    var args = Array.prototype.slice.call(arguments, 1);
    while (i--) arr[length - 1 - i] = createArray.apply(this, args);
  }

  return arr;
}

function pitchBetween(a, b) {
  let x = a.x - b.x;
  let y = a.y - b.y;
  let z = a.z - b.z;
  let r = sqrt(x * x + y * y + z * z);
  let result = atan(z / x);
  if (x > 0) return result;
  else if (x < 0) return result + 180;
  else return 0;
}

function azimuthBetween(a, b) {
  let x = a.x - b.x;
  let y = a.y - b.y;
  let z = a.z - b.z;
  let r = sqrt(x * x + y * y + z * z);
  let result = acos(y / r);
  if (r = !0) return result;
  else return 0;
}

//////////////////////////////////////////////////

class startingConditions {
  constructor(x, y, z, xVelocity, yVelocity, zVelocity) {
    this.x = x;
    this.xVelocity = xVelocity;
    this.y = y;
    this.yVelocity = yVelocity;
    this.z = z;
    this.zVelocity = zVelocity;
  }
}

class Rocket {
  constructor(c, _maxThrust, _wetMass, _dryMass, _fuelFlow) {
    this.maxThrust = _maxThrust;
    this.wetMass = _wetMass;
    this.dryMass = _dryMass;
    this.fuelFlow = _fuelFlow;

    this.x = c.x;
    this.y = c.y;
    this.z = c.z;
    this.xVelocity = c.xVelocity;
    this.yVelocity = c.yVelocity;
    this.zVelocity = c.zVelocity;

    this.mass = this.wetMass;

    this.xAngle = 0;
    this.yAngle = 0;
    this.zAngle = 0;

    this.time = 0;
    this.throttle = 0;

    this.xEngineForce = 0;
    this.yEngineForce = 0;
    this.zEngineForce = 0;

    this.xGravityForce = 0;
    this.yGravityForce = 0;
    this.zGravityForce = 0;

    this.xForce = 0;
    this.yForce = 0;
    this.zForce = 0;

    this.xAcceleration = 0;
    this.yAcceleration = 0;
    this.zAcceleration = 0;

    this.xGravityAngle = 0;
    this.yGravityAngle = 0;
    this.zGravityAngle = 0;

    this.projected;
    this.projection = [
      [1, 0, 0],
      [0, 1, 0]
    ];
  }

  simulate(_interval) {
    if (_interval > 0) {
      //////////    General Variables    \\\\\\\\\\
      this.time += _interval;

      // Angle to gravity
      this.gravityAzimuth = azimuthBetween(new vector(this.x, this.y, this.z), new vector(0, 0, 0));
      this.gravityPitch = pitchBetween(new vector(this.x, this.y, this.z), new vector(0, 0, 0));
      let gravityX = sin(this.gravityAzimuth) * cos(this.gravityPitch);
      let gravityY = cos(this.gravityAzimuth);
      let gravityZ = sin(this.gravityAzimuth) * sin(this.gravityPitch);


      //////////    X Dimension    \\\\\\\\\\
      //if (this.mass <= this.dryMass) this.xEngineForce = 0;
      //else this.xEngineForce = this.maxThrust * this.throttle * sin(this.zAngle);
      let distanceToCenter = sqrt(this.x * this.x + this.y * this.y + this.z * this.z);

      this.xGravityForce = -gravityAtDistance(sun, distanceToCenter) * this.mass * gravityX;

      this.xForce = this.xEngineForce + this.xGravityForce;
      this.xAcceleration = this.xForce / this.mass;

      this.xVelocity += this.xAcceleration * _interval;
      this.x += this.xVelocity * _interval + (0.5 * this.xAcceleration * _interval * _interval);

      //////////    Y Dimension    \\\\\\\\\\
      //if (this.mass <= this.dryMass) this.yEngineForce = 0;
      //else this.yEngineForce = this.maxThrust * this.throttle * cos(this.zAngle);
      this.yGravityForce = -gravityAtDistance(sun, distanceToCenter) * this.mass * gravityY;

      this.yForce = this.yEngineForce + this.yGravityForce;
      this.yAcceleration = this.yForce / this.mass;

      this.yVelocity += this.yAcceleration * _interval;
      this.y += this.yVelocity * _interval + (0.5 * (this.yAcceleration) * _interval * _interval);
      //////////    Z Dimension    \\\\\\\\\\
      //if (this.mass <= this.dryMass) this.zEngineForce = 0;
      //else this.zEngineForce = this.maxThrust * this.throttle * cos(this.zAngle);
      this.zGravityForce = -gravityAtDistance(sun, distanceToCenter) * this.mass * gravityZ;

      this.zForce = this.zEngineForce + this.zGravityForce;
      this.zAcceleration = this.zForce / this.mass;

      this.zVelocity += this.zAcceleration * _interval;
      this.z += this.zVelocity * _interval + (0.5 * (this.zAcceleration) * _interval * _interval);


      //////////    Simulation Handling    \\\\\\\\\\
      if (this.mass <= this.dryMass) this.mass = this.dryMass;
      else this.mass -= (this.fuelFlow * _interval * this.throttle);
    }
  }

  draw() {
    point3d(this.x, this.y, this.z);
    push();
    noStroke();
    text3d(round(this.x), this.x + 10 * zoom, this.y + -20 * zoom, this.z);
    text3d(round(this.y), this.x + 10 * zoom, this.y, this.z);
    text3d(round(this.z), this.x + 10 * zoom, this.y + 20 * zoom, this.z);
    pop();
  }

  project() {
    this.projected = matrixToVector(matrixMultiplyVector(this.projection, this));
    this.projected.z = 0;
  }
}

class Body {
  constructor(x, y, z, mass, radius, name) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.mass = mass;
    this.radius = radius;
    this.name = name;
    this.vector = new vector(x, y, z);
  }

  draw() {
    ellipse3d(this.vector.projected.x, this.vector.projected.y, this.vector.projected.z, this.radius * 2, this.radius * 2);
    text3d(this.name, this.vector.projected.x + this.radius, this.vector.projected.y, this.vector.projected.z);
  }

  drawUnscaled(size) {
    ellipse3d(this.vector.projected.x, this.vector.projected.y, this.vector.projected.z, size * this.radius * 2 * zoom, size * this.radius * 2 * zoom);
    text3d(this.name, this.vector.projected.x + this.radius, this.vector.projected.y, this.vector.projected.z);
  }

  drawAbsolute(size) {
    ellipse3d(this.vector.projected.x, this.vector.projected.y, this.vector.projected.z, size * zoom, size * zoom);
    text3d(this.name, this.vector.projected.x + this.radius, this.vector.projected.y, this.vector.projected.z);
  }

  rotate(xAngle, yAngle, zAngle) {
    this.vector.rotate(xAngle, yAngle, zAngle);
  }

  project() {
    this.vector.project();
  }

  rotateX(angle) {
    this.vector.rotateX(angle);
  }

  rotateY(angle) {
    this.vector.rotateY(angle);
  }

  rotateZ(angle) {
    this.vector.rotateZ(angle);
  }
}

function gravityAtDistance(body, distance) {
  return ((body.mass * gravitationalConstant) / (distance * distance));
}

function airResistance(dragCoefficient, airDensity, velocity, referenceArea) {
  return (0.5 * dragCoefficient * airDensity * velocity * velocity * referenceArea);
}

function airDensity(airPressure, temperature) {
  return (airPressure * 100 / (287.058 * (temperature + 273.15)));
}

function orbitalVelocity(G, mass, radius) {
  return sqrt((G * mass) / radius);
}
