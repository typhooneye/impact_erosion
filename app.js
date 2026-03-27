const CONSTANTS = {
  G: 6.67e-11,
  KB: 1.38e-23,
  MH: 1.674e-27,
  EARTH_MASS: 5.972e24,
  EARTH_RADIUS: 6.378e6,
};

const state = {
  planetesimalMode: "single",
};

const coeffSets = {
  shuvalov: { A12: 1.4, Afr: 1.4, B12: 5.5, Bfr: 3.5, label: "Shuvalov et al. (2014)" },
  sinclair: { A12: 1.4, Afr: 1.4, B12: 1.24, Bfr: 1.97, label: "Sinclair & Wyatt (2021)" },
};

function $(id) {
  return document.getElementById(id);
}

function readNumber(id) {
  const value = Number($(id).value);
  return Number.isFinite(value) ? value : NaN;
}

function clamp(value, min, max) {
  return Math.min(max, Math.max(min, value));
}

function clampPositive(value) {
  if (!Number.isFinite(value)) return 0;
  return Math.max(0, value);
}

function safeLog10(value) {
  return Math.log(value) / Math.log(10);
}

function formatNumber(value, digits = 3) {
  if (!Number.isFinite(value)) return "not defined";
  if (value === 0) return "0";
  const abs = Math.abs(value);
  if (abs >= 1e4 || abs < 1e-3) return value.toExponential(digits);
  return value.toLocaleString(undefined, { maximumFractionDigits: digits });
}

function formatMassKg(value) {
  return `${formatNumber(value, 4)} kg`;
}

function formatMeters(value) {
  if (!Number.isFinite(value)) return "not defined";
  if (Math.abs(value) >= 1000) return `${formatNumber(value / 1000, 4)} km`;
  return `${formatNumber(value, 4)} m`;
}

function buildPlanet() {
  const planetMassEarth = readNumber("planetMassEarth");
  const mu = readNumber("meanMolecularWeight");
  const rhoAtm = readNumber("atmosDensity");
  const temperature = readNumber("atmosTemperature");

  if ([planetMassEarth, mu, rhoAtm, temperature].some((x) => !Number.isFinite(x) || x <= 0)) {
    throw new Error("Planet inputs must all be positive numbers.");
  }

  const mass = planetMassEarth * CONSTANTS.EARTH_MASS;
  const radius = CONSTANTS.EARTH_RADIUS * Math.pow(mass / CONSTANTS.EARTH_MASS, 0.27);
  const gravity = (CONSTANTS.G * mass) / (radius ** 2);
  const scaleHeight = (CONSTANTS.KB * temperature) / (mu * CONSTANTS.MH * gravity);
  const density = mass / ((4 / 3) * Math.PI * radius ** 3);
  const escapeVelocity = Math.sqrt((2 * CONSTANTS.G * mass) / radius);
  const atmosphereMass = ((4 / 3) * Math.PI * ((radius + scaleHeight) ** 3 - radius ** 3)) * rhoAtm;
  const fullVolume = (4 / 3) * Math.PI * radius ** 3;

  return {
    mass,
    radius,
    gravity,
    scaleHeight,
    density,
    escapeVelocity,
    atmosphereMass,
    rhoAtm,
    fullVolume,
  };
}

function impactorDiameterFromMass(mass, density) {
  return Math.cbrt((6 * mass) / (density * Math.PI));
}

function chiFromXi(xi) {
  if (!Number.isFinite(xi) || xi <= 0) return 0;
  if (xi > 1e6) return 10 ** (-0.6438 * xi + 0.4746);

  const logXi = safeLog10(xi);
  const exponent =
    -6.375 +
    5.239 * logXi -
    2.121 * logXi ** 2 +
    0.397 * logXi ** 3 -
    0.037 * logXi ** 4 +
    0.0013 * logXi ** 5;
  return 10 ** exponent;
}

function interpolateFragmentationLoss(crateringLoss, burstLoss, H12, Hfr) {
  if (crateringLoss <= 0 || burstLoss <= 0 || Hfr === H12) return Math.max(crateringLoss, burstLoss, 0);
  const weight = Hfr / (Hfr - H12);
  const logLoss = Math.log(crateringLoss) + (Math.log(burstLoss) - Math.log(crateringLoss)) * weight;
  return clampPositive(Math.exp(logLoss));
}

function computeSinglePlanetesimal(planet, options) {
  const mass = options.mass;
  const density = options.density;
  const velocityRatio = options.velocityRatio;
  const coeff = coeffSets[options.coeffSet];

  if (!Number.isFinite(mass) || mass <= 0) throw new Error("Single impactor mass must be positive.");
  if (!Number.isFinite(density) || density <= 0) throw new Error("Impactor density must be positive.");
  if (!Number.isFinite(velocityRatio) || velocityRatio < 0) throw new Error("Velocity ratio must be non-negative.");

  const diameter = impactorDiameterFromMass(mass, density);
  const velocity = velocityRatio * planet.escapeVelocity;
  const xi =
    (diameter / planet.scaleHeight) ** 3 *
    (velocityRatio ** 2 - 1) *
    ((density * planet.density) / (planet.rhoAtm * (density + planet.density)));

  const chi = chiFromXi(xi);
  const crateringCap = 2 * Math.PI * planet.rhoAtm * planet.scaleHeight ** 2 * planet.radius;
  const crateringLoss = clampPositive(Math.min(crateringCap, chi * (velocityRatio ** 2 - 1) * mass));

  const densityRatio = Math.sqrt(density / planet.rhoAtm);
  const H12 = -coeff.A12 * planet.scaleHeight * Math.log((coeff.B12 * diameter / (2 * planet.scaleHeight)) * densityRatio);
  const Hfr = -coeff.Afr * planet.scaleHeight * Math.log((coeff.Bfr * diameter / (2 * planet.scaleHeight)) * densityRatio);

  const rhoImpGcm = density / 1000;
  const velocityKm = velocity / 1000;
  const scaleHeightKm = planet.scaleHeight / 1000;
  const diameterKm = diameter / 1000;
  const burstLoss = clampPositive(
    5.5e-3 * (velocityKm - 15 - 0.2 * scaleHeightKm) ** 2 * Math.sqrt(rhoImpGcm) * (diameterKm / scaleHeightKm) * mass
  );

  let regime = "transitional";
  let atmosphericLoss = Math.max(crateringLoss, burstLoss);

  if (H12 < 0 && Hfr < 0) {
    regime = "cratering";
    atmosphericLoss = crateringLoss;
  } else if (H12 > 0 && Hfr > 0) {
    regime = "aerial burst";
    atmosphericLoss = burstLoss;
  } else if (H12 < 0 && Hfr > 0) {
    regime = "fragmentation";
    atmosphericLoss = interpolateFragmentationLoss(crateringLoss, burstLoss, H12, Hfr);
  }

  return {
    mass,
    diameter,
    velocity,
    velocityRatio,
    xi,
    chi,
    H12,
    Hfr,
    regime,
    atmosphericLoss,
    lossFraction: clampPositive(atmosphericLoss / planet.atmosphereMass),
    coefficientLabel: coeff.label,
  };
}

function integratePowerLawMass(dMin, dMax, q) {
  const exponent = 4 - q;
  if (Math.abs(exponent) < 1e-12) return Math.log(dMax / dMin);
  return (dMax ** exponent - dMin ** exponent) / exponent;
}

function integrateCountBin(k, left, right, q) {
  const exponent = 1 - q;
  if (Math.abs(exponent) < 1e-12) return k * Math.log(right / left);
  return (k * (right ** exponent - left ** exponent)) / exponent;
}

function logspace(min, max, count) {
  const start = Math.log(min);
  const end = Math.log(max);
  const step = (end - start) / (count - 1);
  const values = [];
  for (let i = 0; i < count; i += 1) values.push(Math.exp(start + step * i));
  return values;
}

function computeCollectivePlanetesimal(planet, options) {
  const totalMass = options.totalMass;
  const density = options.density;
  const velocityRatio = options.velocityRatio;
  const dMin = options.dMin;
  const dMax = options.dMax;
  const q = options.q;
  const bins = options.bins;

  if (!Number.isFinite(totalMass) || totalMass <= 0) throw new Error("Collective mass must be positive.");
  if (!Number.isFinite(density) || density <= 0) throw new Error("Impactor density must be positive.");
  if (!Number.isFinite(velocityRatio) || velocityRatio < 0) throw new Error("Velocity ratio must be non-negative.");
  if (!Number.isFinite(dMin) || !Number.isFinite(dMax) || dMin <= 0 || dMax <= dMin) {
    throw new Error("Diameter bounds must be positive and Dmax must be larger than Dmin.");
  }
  if (!Number.isFinite(q) || q <= 0) throw new Error("Power-law index q must be positive.");
  if (!Number.isFinite(bins) || bins < 20) throw new Error("Use at least 20 bins for the collective integration.");

  const edges = logspace(dMin, dMax, Math.floor(bins) + 1);
  const massIntegral = integratePowerLawMass(dMin, dMax, q);
  const k = totalMass / (((Math.PI / 6) * density) * massIntegral);

  let expectedImpactors = 0;
  let totalLoss = 0;
  let crateringLoss = 0;
  let burstLoss = 0;
  let fragmentationLoss = 0;

  for (let i = 0; i < edges.length - 1; i += 1) {
    const left = edges[i];
    const right = edges[i + 1];
    const mid = Math.sqrt(left * right);
    const count = integrateCountBin(k, left, right, q);
    const mass = (Math.PI / 6) * density * mid ** 3;
    const single = computeSinglePlanetesimal(planet, {
      mass,
      density,
      velocityRatio,
      coeffSet: options.coeffSet,
    });

    expectedImpactors += count;
    totalLoss += count * single.atmosphericLoss;

    if (single.regime === "cratering") crateringLoss += count * single.atmosphericLoss;
    else if (single.regime === "aerial burst") burstLoss += count * single.atmosphericLoss;
    else if (single.regime === "fragmentation") fragmentationLoss += count * single.atmosphericLoss;
  }

  return {
    expectedImpactors,
    totalLoss,
    lossFraction: clampPositive(totalLoss / planet.atmosphereMass),
    crateringLoss,
    burstLoss,
    fragmentationLoss,
    k,
    q,
    dMin,
    dMax,
  };
}

function sphericalCapVolume(radius, height) {
  const cappedHeight = clamp(height, 0, 2 * radius);
  return (Math.PI * cappedHeight ** 2 * (3 * radius - cappedHeight)) / 3;
}

function computeGiantImpact(planet, options) {
  const massEarth = options.massEarth;
  const velocityRatio = options.velocityRatio;
  const impactParameter = options.impactParameter;

  if (!Number.isFinite(massEarth) || massEarth <= 0) throw new Error("Impactor mass must be positive.");
  if (!Number.isFinite(velocityRatio) || velocityRatio < 0) throw new Error("Velocity ratio must be non-negative.");
  if (!Number.isFinite(impactParameter) || impactParameter < 0 || impactParameter > 1) {
    throw new Error("Impact parameter b must be between 0 and 1.");
  }

  const mass = massEarth * CONSTANTS.EARTH_MASS;
  const radius = CONSTANTS.EARTH_RADIUS * Math.pow(mass / CONSTANTS.EARTH_MASS, 0.27);
  const density = mass / ((4 / 3) * Math.PI * radius ** 3);
  const fullVolume = (4 / 3) * Math.PI * radius ** 3;
  const velocity = velocityRatio * planet.escapeVelocity;
  const dCap = Math.max(0, (planet.radius + radius) * (1 - impactParameter));
  const vCapImpactor = sphericalCapVolume(radius, dCap);
  const vCapPlanet = sphericalCapVolume(planet.radius, dCap);
  const fM = (planet.density * vCapPlanet + density * vCapImpactor) / (planet.density * planet.fullVolume + density * fullVolume);
  const ratio =
    0.64 *
    (
      (velocity / planet.escapeVelocity) ** 2 *
      Math.sqrt(mass / (mass + planet.mass)) *
      Math.sqrt(density / planet.density) *
      fM
    ) ** 0.65;

  const erosionFraction = clamp(ratio, 0, 1);
  const atmosphereLoss = erosionFraction * planet.atmosphereMass;

  return {
    radius,
    density,
    velocity,
    impactParameter,
    fM,
    erosionFraction,
    atmosphereLoss,
  };
}

function renderPlanetSummary(planet) {
  $("planetSummary").innerHTML = `
    <strong>Derived planet properties</strong><br />
    Radius: ${formatMeters(planet.radius)}<br />
    Escape velocity: ${formatNumber(planet.escapeVelocity / 1000, 4)} km/s<br />
    Surface gravity: ${formatNumber(planet.gravity, 4)} m/s²<br />
    Scale height: ${formatMeters(planet.scaleHeight)}<br />
    Atmospheric mass estimate: ${formatMassKg(planet.atmosphereMass)}
  `;
}

function resultCard(title, value, subvalue) {
  return `
    <article class="card">
      <h3>${title}</h3>
      <div class="value">${value}</div>
      <div class="subvalue">${subvalue}</div>
    </article>
  `;
}

function renderPlanetesimalResults(data, mode) {
  const results = [];

  if (mode === "single") {
    results.push(resultCard("Regime", data.regime, `Coefficient set: ${data.coefficientLabel}`));
    results.push(resultCard("Atmosphere lost", formatMassKg(data.atmosphericLoss), `${formatNumber(data.lossFraction * 100, 4)}% of estimated atmosphere`));
    results.push(resultCard("Impactor diameter", formatMeters(data.diameter), `Velocity: ${formatNumber(data.velocity / 1000, 4)} km/s`));
    results.push(resultCard("Erosion efficiency ξ", formatNumber(data.xi, 5), `χ(ξ): ${formatNumber(data.chi, 5)}`));
    results.push(resultCard("H1/2", formatMeters(data.H12), "Half-velocity altitude"));
    results.push(resultCard("Hfr", formatMeters(data.Hfr), "Fragmentation altitude"));
  } else {
    results.push(resultCard("Total atmosphere lost", formatMassKg(data.totalLoss), `${formatNumber(data.lossFraction * 100, 4)}% of estimated atmosphere`));
    results.push(resultCard("Expected number of impactors", formatNumber(data.expectedImpactors, 4), `q = ${formatNumber(data.q, 3)}`));
    results.push(resultCard("Cratering contribution", formatMassKg(data.crateringLoss), "Integrated over the full size distribution"));
    results.push(resultCard("Burst contribution", formatMassKg(data.burstLoss), "Integrated over the full size distribution"));
    results.push(resultCard("Fragmentation contribution", formatMassKg(data.fragmentationLoss), "Integrated over the full size distribution"));
    results.push(resultCard("Distribution constant k", formatNumber(data.k, 4), `Diameter range: ${formatMeters(data.dMin)} to ${formatMeters(data.dMax)}`));
  }

  $("planetesimalResults").innerHTML = results.join("");
}

function renderGiantResults(data) {
  $("giantResults").innerHTML = [
    resultCard("Atmosphere lost", formatMassKg(data.atmosphereLoss), `${formatNumber(data.erosionFraction * 100, 4)}% of estimated atmosphere`),
    resultCard("Interacting mass fraction fM", formatNumber(data.fM, 5), `Impact parameter b = ${formatNumber(data.impactParameter, 3)}`),
    resultCard("Impactor radius", formatMeters(data.radius), `Impactor density: ${formatNumber(data.density, 4)} kg/m³`),
    resultCard("Impact velocity", `${formatNumber(data.velocity / 1000, 4)} km/s`, "Scaled from the target escape velocity"),
  ].join("");
}

function showError(targetId, message) {
  $(targetId).innerHTML = resultCard("Input error", "Check values", message);
}

function updatePlanetesimalMode(mode) {
  state.planetesimalMode = mode;
  $("singleInputs").classList.toggle("is-hidden", mode !== "single");
  $("collectiveInputs").classList.toggle("is-hidden", mode !== "collective");
  $("singleModeBtn").classList.toggle("is-active", mode === "single");
  $("collectiveModeBtn").classList.toggle("is-active", mode === "collective");
}

function calculatePlanetesimal() {
  try {
    const planet = buildPlanet();
    renderPlanetSummary(planet);

    const shared = {
      density: readNumber("planetesimalDensity"),
      velocityRatio: readNumber("planetesimalVelocityRatio"),
      coeffSet: $("planetesimalCoeffSet").value,
    };

    if (state.planetesimalMode === "single") {
      const result = computeSinglePlanetesimal(planet, {
        ...shared,
        mass: readNumber("singleImpactorMass"),
      });
      renderPlanetesimalResults(result, "single");
    } else {
      const result = computeCollectivePlanetesimal(planet, {
        ...shared,
        totalMass: readNumber("collectiveMass"),
        dMin: readNumber("diameterMin"),
        dMax: readNumber("diameterMax"),
        q: readNumber("powerLawIndex"),
        bins: readNumber("collectiveBins"),
      });
      renderPlanetesimalResults(result, "collective");
    }
  } catch (error) {
    showError("planetesimalResults", error.message);
  }
}

function calculateGiant() {
  try {
    const planet = buildPlanet();
    renderPlanetSummary(planet);

    const result = computeGiantImpact(planet, {
      massEarth: readNumber("giantMassEarth"),
      velocityRatio: readNumber("giantVelocityRatio"),
      impactParameter: readNumber("impactParameter"),
    });
    renderGiantResults(result);
  } catch (error) {
    showError("giantResults", error.message);
  }
}

function bindEvents() {
  $("singleModeBtn").addEventListener("click", () => updatePlanetesimalMode("single"));
  $("collectiveModeBtn").addEventListener("click", () => updatePlanetesimalMode("collective"));
  $("calculatePlanetesimal").addEventListener("click", calculatePlanetesimal);
  $("calculateGiant").addEventListener("click", calculateGiant);

  ["planetMassEarth", "meanMolecularWeight", "atmosDensity", "atmosTemperature"].forEach((id) => {
    $(id).addEventListener("input", () => {
      try {
        renderPlanetSummary(buildPlanet());
      } catch {
        $("planetSummary").textContent = "Enter valid planet inputs to see derived quantities.";
      }
    });
  });
}

function init() {
  bindEvents();
  renderPlanetSummary(buildPlanet());
  calculatePlanetesimal();
  calculateGiant();
}

init();
