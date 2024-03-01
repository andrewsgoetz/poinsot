import * as THREE from 'three';
import { GUI } from 'lil-gui';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const ZERO_VECTOR = new THREE.Vector3(0., 0., 0.);
const arrowHeadRelativeScale = 0.03;

// Initial values for user-controllable inputs.
const initialMass = 100.;
const initialA = 0.1;
const initialB = 0.2;
const initialC = 0.3;
const initialAngularMomentum = 3.;
const initialQ0unnormalized = 1.;
const initialQ1unnormalized = -0.7;
const initialQ2unnormalized = 0.4;
const initialQ3unnormalized = -1.;

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

const axesHelper = new THREE.AxesHelper();

const omegaArrow = new THREE.ArrowHelper(new THREE.Vector3(1., 0., 0.), ZERO_VECTOR, 1., 0x00ffff);

const angularMomentumArrow = new THREE.ArrowHelper(
    new THREE.Vector3(0., 0., -1.),
    ZERO_VECTOR,
    initialAngularMomentum,
    0xffffff,
    arrowHeadRelativeScale * initialAngularMomentum,
    arrowHeadRelativeScale * initialAngularMomentum,
);

const rigidBodyGeometry = new THREE.SphereGeometry(1.);
const rigidBodyWireframe = new THREE.WireframeGeometry(rigidBodyGeometry);
const rigidBodyObject = new THREE.LineSegments(rigidBodyWireframe);
rigidBodyObject.material.color.set(0xffffff);
rigidBodyObject.scale.set(initialA, initialB, initialC);

const momentalEllipsoidGeometry = new THREE.SphereGeometry(1.);
const momentalEllipsoidWireframe = new THREE.WireframeGeometry(momentalEllipsoidGeometry);
const momentalEllipsoidMaterial = new THREE.LineBasicMaterial({ color: 0x666666 });
const momentalEllipsoidObject = new THREE.LineSegments(momentalEllipsoidWireframe, momentalEllipsoidMaterial);

const invariablePlane = new THREE.Plane(new THREE.Vector3(0., 0., 1.), 1.);
const invariablePlaneHelper = new THREE.PlaneHelper(invariablePlane, 10., 0xffff00);

let polhodePoints;
const polhodeGeometry = new THREE.BufferGeometry();
const polhodeMaterial = new THREE.LineBasicMaterial({ color: 0xff00ff, linewidth: 2, });
const polhodeObject = new THREE.Line(polhodeGeometry, polhodeMaterial);

let herpolhodePoints;
const herpolhodeGeometry = new THREE.BufferGeometry();
const herpolhodeMaterial = new THREE.LineBasicMaterial({ color: 0xffffff, linewidth: 2, });
const herpolhodeObject = new THREE.Line(herpolhodeGeometry, herpolhodeMaterial);

scene.add(
    axesHelper,
    angularMomentumArrow,
    omegaArrow,
    rigidBodyObject,
    momentalEllipsoidObject,
    invariablePlaneHelper,
    polhodeObject,
    herpolhodeObject,
);

let inputsFrozen = false;
let animating = false;

let timeStampOfLastPointsSaved;
let previousTimeStamp;
function animate(timeStamp) {
    if (animating) {
        if (!previousTimeStamp) {
            previousTimeStamp = timeStamp;
        } else {
            const I1 = ellipsoidProperties.I1;
            const I2 = ellipsoidProperties.I2;
            const I3 = ellipsoidProperties.I3;
            let q0 = quaternion.q0;
            let q1 = quaternion.q1;
            let q2 = quaternion.q2;
            let q3 = quaternion.q3;
            let omega1b = angularVelocity.omega1b;
            let omega2b = angularVelocity.omega2b;
            let omega3b = angularVelocity.omega3b;

            const result = integrate(I1, I2, I3, q0, q1, q2, q3, omega1b, omega2b, omega3b, previousTimeStamp, timeStamp);
            previousTimeStamp = timeStamp;

            quaternion.q0 = result.q0;
            quaternion.q1 = result.q1;
            quaternion.q2 = result.q2;
            quaternion.q3 = result.q3;
            quaternion.q = new THREE.Quaternion(q1, q2, q3, q0);
            updateQuaternion();

            angularVelocity.omega1b = result.omega1b;
            angularVelocity.omega2b = result.omega2b;
            angularVelocity.omega3b = result.omega3b;
            const omega = new THREE.Vector3(omega1b, omega2b, omega3b).applyQuaternion(quaternion.q);
            angularVelocity.omega1 = omega.x;
            angularVelocity.omega2 = omega.y;
            angularVelocity.omega3 = omega.z;

            rigidBodyObject.setRotationFromQuaternion(quaternion.q);
            momentalEllipsoidObject.setRotationFromQuaternion(quaternion.q);
            updateAngularVelocityVector();

            polhodeGeometry.setFromPoints(polhodePoints);
            herpolhodeGeometry.setFromPoints(herpolhodePoints);
            polhodeObject.setRotationFromQuaternion(quaternion.q);
        }
    }
    renderer.render(scene, camera);
    requestAnimationFrame(animate);
}
window.requestAnimationFrame(animate);

const animation = {
    toggle: toggleAnimation,
    reset: resetAnimation,
    info: showOrHideInfo,
}

const ellipsoidProperties = {
    M: initialMass,
    a: initialA,
    b: initialB,
    c: initialC,
    I1: 0.2 * initialMass * (initialB * initialB + initialC * initialC),
    I2: 0.2 * initialMass * (initialC * initialC + initialA * initialA),
    I3: 0.2 * initialMass * (initialA * initialA + initialB * initialB),
};

const quaternion = {
    q0unnormalized: initialQ0unnormalized,
    q1unnormalized: initialQ1unnormalized,
    q2unnormalized: initialQ2unnormalized,
    q3unnormalized: initialQ3unnormalized,
    // values below set during initialization
    q0: 1.,
    q1: 0.,
    q2: 0.,
    q3: 0.,
    q: new THREE.Quaternion(0., 0., 0., 1.),
}

// all values set during initialization
const angularVelocity = {
    // inertial frame
    omega1: 1.,
    omega2: 0.,
    omega3: 0.,
    // body frame
    omega1b: 1.,
    omega2b: 0.,
    omega3b: 0.,
}

const constantsOfMotion = {
    kineticEnergy: 0., // set during initialization
    angularMomentum: initialAngularMomentum,
}

const graphicsProperties = {
    showObject: true,
    showAxes: true,
    showMomentalEllipsoid: true,
    showAngularMomentumVector: true,
    showAngularVelocityVector: true,
    showInvariablePlane: true,
    showPolhode: true,
    showHerpolhode: true,
}

const gui = new GUI();

const ellipsoidFolder = gui.addFolder('Ellipsoid Rigid Body');
const massController = ellipsoidFolder.add(ellipsoidProperties, 'M', 0, 500).name('M').onChange(() => userUpdatedEllipsoid());
const aController = ellipsoidFolder.add(ellipsoidProperties, 'a', 0, 1).name('a').onChange(() => userUpdatedEllipsoid());
const bController = ellipsoidFolder.add(ellipsoidProperties, 'b', 0, 1).name('b').onChange(() => userUpdatedEllipsoid());
const cController = ellipsoidFolder.add(ellipsoidProperties, 'c', 0, 1).name('c').onChange(() => userUpdatedEllipsoid());

const quaternionFolder = gui.addFolder('Quaternion Orientation');
const q0Controller = quaternionFolder.add(quaternion, 'q0unnormalized', -1., 1.).name('q<sub>0</sub>').onChange(() => userUpdatedQuaternion());
const q1Controller = quaternionFolder.add(quaternion, 'q1unnormalized', -1., 1.).name('q<sub>1</sub>').onChange(() => userUpdatedQuaternion());
const q2Controller = quaternionFolder.add(quaternion, 'q2unnormalized', -1., 1.).name('q<sub>2</sub>').onChange(() => userUpdatedQuaternion());
const q3Controller = quaternionFolder.add(quaternion, 'q3unnormalized', -1., 1.).name('q<sub>3</sub>').onChange(() => userUpdatedQuaternion());

const angularVelocityFolder = gui.addFolder('Angular Velocity');
const omega1Controller = angularVelocityFolder.add(angularVelocity, 'omega1').name('\u03C9<sub>1</sub>').disable();
const omega2Controller = angularVelocityFolder.add(angularVelocity, 'omega2').name('\u03C9<sub>2</sub>').disable();
const omega3Controller = angularVelocityFolder.add(angularVelocity, 'omega3').name('\u03C9<sub>3</sub>').disable();

const constantsOfMotionFolder = gui.addFolder('Constants of Motion');
const angularMomentumController = constantsOfMotionFolder.add(constantsOfMotion, 'angularMomentum', 0, 10).name('L').onChange(() => userUpdatedAngularMomentum());
const kineticEnergyController = constantsOfMotionFolder.add(constantsOfMotion, 'kineticEnergy').name('T').disable();

const toggleGraphicsOptionsFolder = gui.addFolder('Show');
toggleGraphicsOptionsFolder.add(graphicsProperties, 'showObject').name('Object').onChange(v => { rigidBodyObject.visible = v; });
toggleGraphicsOptionsFolder.add(graphicsProperties, 'showAxes').name('Axes').onChange(v => { axesHelper.visible = v; });
toggleGraphicsOptionsFolder.add(graphicsProperties, 'showMomentalEllipsoid').name('Momental Ellipsoid').onChange(v => { momentalEllipsoidObject.visible = v; });
toggleGraphicsOptionsFolder.add(graphicsProperties, 'showAngularMomentumVector').name('Angular Momentum').onChange(v => { angularMomentumArrow.visible = v; });
toggleGraphicsOptionsFolder.add(graphicsProperties, 'showAngularVelocityVector').name('Angular Velocity').onChange(v => { omegaArrow.visible = v; });
toggleGraphicsOptionsFolder.add(graphicsProperties, 'showInvariablePlane').name('Invariable Plane').onChange(v => { invariablePlaneHelper.visible = v; });
toggleGraphicsOptionsFolder.add(graphicsProperties, 'showPolhode').name('Polhode').onChange(v => { polhodeObject.visible = v; });
toggleGraphicsOptionsFolder.add(graphicsProperties, 'showHerpolhode').name('Herpolhode').onChange(v => { herpolhodeObject.visible = v; });

const toggleAnimationButton = gui.add(animation, 'toggle').name(''); // name set elsewhere
const resetAnimationButton = gui.add(animation, 'reset').name('Reset');
const infoButton = gui.add(animation, 'info').name('Hide Info');

const userInputs = [
    massController,
    aController,
    bController,
    cController,
    q0Controller,
    q1Controller,
    q2Controller,
    q3Controller,
    angularMomentumController,
];

init();
const orbitControls = new OrbitControls(camera, renderer.domElement);
orbitControls.enablePan = false;
toggleAnimation();

// Initialize.
function init() {
    userUpdatedEllipsoid();
    userUpdatedQuaternion();
    userUpdatedAngularVelocity();
    userUpdatedAngularMomentum();
    polhodePoints = [];
    herpolhodePoints = [];
    polhodeGeometry.setFromPoints(polhodePoints);
    herpolhodeGeometry.setFromPoints(herpolhodePoints);
    polhodeObject.setRotationFromQuaternion(quaternion.q);
    camera.position.set(-1., 2.5, 2.);
    camera.up = new THREE.Vector3(0., 0., 1.);
    camera.lookAt(ZERO_VECTOR);
    for (const userInput of userInputs) {
        userInput.enable();
        userInput.updateDisplay();
    }
    for (const child of toggleGraphicsOptionsFolder.children) {
        child.setValue(true); // make all objects visible
    }
}

function toggleAnimation() {
    if (!inputsFrozen) {
        inputsFrozen = true;
        for (const userInput of userInputs) {
            userInput.disable();
        }
    }
    if (animating) {
        animating = false;
        previousTimeStamp = undefined;
        toggleAnimationButton.name('Start Animation');
        resetAnimationButton.enable();
    } else {
        animating = true;
        toggleAnimationButton.name('Stop Animation');
        resetAnimationButton.disable();
    }
}

function resetAnimation() {
    ellipsoidProperties.M = initialMass;
    ellipsoidProperties.a = initialA;
    ellipsoidProperties.b = initialB;
    ellipsoidProperties.c = initialC;
    quaternion.q0unnormalized = initialQ0unnormalized;
    quaternion.q1unnormalized = initialQ1unnormalized;
    quaternion.q2unnormalized = initialQ2unnormalized;
    quaternion.q3unnormalized = initialQ3unnormalized;
    constantsOfMotion.angularMomentum = initialAngularMomentum;
    inputsFrozen = false;
    init();
}

function showOrHideInfo() {
    const info = document.getElementById('info');
    console.log(info.style.display);
    if (info.style.display === 'none') {
        info.style.display = 'block';
        infoButton.name('Hide Info');
    } else {
        info.style.display = 'none';
        infoButton.name('Show Info');
    }
}

function userUpdatedEllipsoid() {
    const M = ellipsoidProperties.M;
    const a = ellipsoidProperties.a;
    const b = ellipsoidProperties.b;
    const c = ellipsoidProperties.c;
    ellipsoidProperties.I1 = 0.2 * M * (b * b + c * c);
    ellipsoidProperties.I2 = 0.2 * M * (c * c + a * a);
    ellipsoidProperties.I3 = 0.2 * M * (a * a + b * b);
    rigidBodyObject.scale.set(a, b, c);
    userUpdatedAngularVelocity();
}

function userUpdatedQuaternion() {
    const q0u = quaternion.q0unnormalized;
    const q1u = quaternion.q1unnormalized;
    const q2u = quaternion.q2unnormalized;
    const q3u = quaternion.q3unnormalized;
    const norm = Math.sqrt(q0u * q0u + q1u * q1u + q2u * q2u + q3u * q3u);

    let q0;
    let q1;
    let q2;
    let q3;
    if (norm === 0.) {
        q0 = 1.;
        q1 = 0.;
        q2 = 0.;
        q3 = 0.;
    } else {
        q0 = q0u / norm;
        q1 = q1u / norm;
        q2 = q2u / norm;
        q3 = q3u / norm;
    }

    quaternion.q0 = q0;
    quaternion.q1 = q1;
    quaternion.q2 = q2;
    quaternion.q3 = q3;
    quaternion.q = new THREE.Quaternion(q1, q2, q3, q0);

    rigidBodyObject.setRotationFromQuaternion(quaternion.q);
    momentalEllipsoidObject.setRotationFromQuaternion(quaternion.q);
    userUpdatedAngularVelocity();
}

function userUpdatedAngularMomentum() {
    const L = constantsOfMotion.angularMomentum;
    angularMomentumArrow.setLength(L, arrowHeadRelativeScale * L, arrowHeadRelativeScale * L);
    userUpdatedAngularVelocity();
}

function userUpdatedAngularVelocity() { // This only happens indirectly.
    const L = constantsOfMotion.angularMomentum;
    const Lv = new THREE.Vector3(0., 0., -L);
    const Lb = Lv.applyQuaternion(quaternion.q.clone().conjugate());
    const L1b = Lb.x;
    const L2b = Lb.y;
    const L3b = Lb.z;
    const omega1b = L1b / ellipsoidProperties.I1;
    const omega2b = L2b / ellipsoidProperties.I2;
    const omega3b = L3b / ellipsoidProperties.I3;
    angularVelocity.omega1b = omega1b;
    angularVelocity.omega2b = omega2b;
    angularVelocity.omega3b = omega3b;
    const omega = new THREE.Vector3(omega1b, omega2b, omega3b).applyQuaternion(quaternion.q);
    angularVelocity.omega1 = omega.x;
    angularVelocity.omega2 = omega.y;
    angularVelocity.omega3 = omega.z;
    updateAngularVelocityVector();
    userUpdatedKineticEnergy();
}

function userUpdatedKineticEnergy() { // This only happens indirectly.
    const I1 = ellipsoidProperties.I1;
    const I2 = ellipsoidProperties.I2;
    const I3 = ellipsoidProperties.I3;
    const omega1b = angularVelocity.omega1b;
    const omega2b = angularVelocity.omega2b;
    const omega3b = angularVelocity.omega3b;
    const T = 0.5 * (I1 * omega1b * omega1b + I2 * omega2b * omega2b + I3 * omega3b * omega3b);
    kineticEnergyController.setValue(T);
    userUpdatedInvariablePlane();
    userUpdatedMomentalEllipsoid();
}

function userUpdatedInvariablePlane() { // This only happens indirectly.
    const L = constantsOfMotion.angularMomentum;
    const T = constantsOfMotion.kineticEnergy;
    const h = 2 * T / L;
    invariablePlane.constant = h;
}

function userUpdatedMomentalEllipsoid() { // This only happens indirectly.
    const I1 = ellipsoidProperties.I1;
    const I2 = ellipsoidProperties.I2;
    const I3 = ellipsoidProperties.I3;
    const T = constantsOfMotion.kineticEnergy;
    const a = Math.sqrt(2. * T / I1);
    const b = Math.sqrt(2. * T / I2);
    const c = Math.sqrt(2. * T / I3);
    momentalEllipsoidObject.scale.set(a, b, c);
}

function updateQuaternion() {
    quaternion.q0unnormalized = quaternion.q0;
    quaternion.q1unnormalized = quaternion.q1;
    quaternion.q2unnormalized = quaternion.q2;
    quaternion.q3unnormalized = quaternion.q3;
    q0Controller.updateDisplay();
    q1Controller.updateDisplay();
    q2Controller.updateDisplay();
    q3Controller.updateDisplay();
}

function updateAngularVelocityVector() {
    const omega1 = angularVelocity.omega1;
    const omega2 = angularVelocity.omega2;
    const omega3 = angularVelocity.omega3;
    omega1Controller.updateDisplay();
    omega2Controller.updateDisplay();
    omega3Controller.updateDisplay();
    const length = Math.sqrt(omega1 * omega1 + omega2 * omega2 + omega3 * omega3);
    const direction = new THREE.Vector3(omega1, omega2, omega3).multiplyScalar(1. / length);
    omegaArrow.setDirection(direction);
    omegaArrow.setLength(length, arrowHeadRelativeScale * length, arrowHeadRelativeScale * length);
}

function integrate(I1, I2, I3, q0, q1, q2, q3, omega1b, omega2b, omega3b, previousTimeStamp, timeStamp) {
    const deltaTimeStamp = (timeStamp - previousTimeStamp) / 1000.; // seconds
    const maxDt = 0.001; // seconds
    const n = Math.ceil(deltaTimeStamp / maxDt);
    const dt = deltaTimeStamp / n;
    const dtMillis = 1000. * dt;
    for (let i = 0; i < n; ++i) {
        const q0Dot = +0.5 * (-q1 * omega1b - q2 * omega2b - q3 * omega3b);
        const q1Dot = +0.5 * (+q0 * omega1b - q3 * omega2b + q2 * omega3b);
        const q2Dot = +0.5 * (+q3 * omega1b + q0 * omega2b - q1 * omega3b);
        const q3Dot = +0.5 * (-q2 * omega1b + q1 * omega2b + q0 * omega3b);
        q0 += q0Dot * dt;
        q1 += q1Dot * dt;
        q2 += q2Dot * dt;
        q3 += q3Dot * dt;

        // Keep quaternion normalized.
        const qNorm = Math.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
        q0 /= qNorm;
        q1 /= qNorm;
        q2 /= qNorm;
        q3 /= qNorm;

        const omega1bDot = ((I2 - I3) / I1) * omega2b * omega3b;
        const omega2bDot = ((I3 - I1) / I2) * omega3b * omega1b;
        const omega3bDot = ((I1 - I2) / I3) * omega1b * omega2b;
        omega1b += omega1bDot * dt;
        omega2b += omega2bDot * dt;
        omega3b += omega3bDot * dt;

        // Keep kinetic energy constant.
        const T = constantsOfMotion.kineticEnergy;
        const T2 = 0.5 * (I1 * omega1b * omega1b + I2 * omega2b * omega2b + I3 * omega3b * omega3b);
        const r = Math.sqrt(T / T2);
        omega1b *= r;
        omega2b *= r;
        omega3b *= r;

        const currentTimeMillis = previousTimeStamp + i * dtMillis;
        if (!timeStampOfLastPointsSaved || currentTimeMillis - timeStampOfLastPointsSaved > 50) {
            const omega = new THREE.Vector3(omega1b, omega2b, omega3b).applyQuaternion(new THREE.Quaternion(q1, q2, q3, q0));
            const omega1 = omega.x;
            const omega2 = omega.y;
            const omega3 = omega.z;
            polhodePoints.push(new THREE.Vector3(omega1b, omega2b, omega3b));
            herpolhodePoints.push(new THREE.Vector3(omega1, omega2, omega3));
            timeStampOfLastPointsSaved = currentTimeMillis;
        }
    }
    return { q0, q1, q2, q3, omega1b, omega2b, omega3b };
}
