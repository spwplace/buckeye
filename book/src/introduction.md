# Introduction

> "The release of atom power has changed everything except our way of thinking."
> — Albert Einstein

## The Question

The "nuclear football" is one of the most famous objects in American politics—a 45-pound aluminum briefcase officially called the "Presidential Emergency Satchel," carried at all times by a military aide within arm's reach of the President. It contains the communication equipment and authentication codes needed to authorize a nuclear strike.

But what if we took the term literally?

**Can an actual NFL football, filled with weapons-grade plutonium, achieve nuclear criticality?**

This is not a frivolous question. Answering it requires us to understand:

1. The physics of nuclear fission and neutron chain reactions
2. The properties of fissile materials like plutonium-239
3. The mathematics of neutron transport
4. The geometry of critical assemblies
5. Monte Carlo computational methods

This book will take you through all of it.

## Why This Matters

Beyond the admittedly absurd premise, this project is a serious exercise in computational nuclear physics. The same Monte Carlo methods we use here are employed by national laboratories around the world to:

- Design nuclear reactors
- Assess weapons performance
- Ensure criticality safety in fuel processing facilities
- Plan radiation shielding
- Develop medical isotope production

Understanding these methods gives insight into some of the most carefully guarded and consequential physics of the modern era.

## The Answer (Spoiler)

**Yes.** An NFL football filled with delta-phase plutonium-239 would be a supercritical nuclear assembly.

Our Monte Carlo simulation calculates k_eff ≈ 1.2 for a Pu-filled football, meaning each generation of neutrons produces about 20% more neutrons than the previous generation. The neutron population would grow exponentially—this is a prompt critical excursion.

The football shape is suboptimal (a sphere minimizes neutron leakage), but the ~22 kg of plutonium in an NFL-sized football is roughly twice the critical mass. Even with the pointed-end geometry causing extra leakage, there's plenty of margin for criticality.

## What We'll Cover

### Part I: Nuclear Physics Fundamentals
We start with the nucleus itself—what holds it together, why some nuclei are stable and others fission, and what happens when a neutron strikes plutonium-239.

### Part II: Monte Carlo Neutron Transport
The heart of our simulation. We'll derive the algorithms used to track individual neutrons through matter, including how to sample random collision distances, reaction types, and scattering angles.

### Part III: Cross Sections
The probability of nuclear reactions, and how it varies with neutron energy. We'll cover the ENDF data format and how to interpolate continuous-energy cross sections.

### Part IV: Geometry
Ray tracing through complex shapes. We'll derive the mathematics for spheres, ellipsoids, and the superellipsoid that models a football's pointed ends.

### Part V: Implementation
A walkthrough of the Julia code, connecting the physics and mathematics to actual implementation.

### Part VI: Results
Validation against the Jezebel benchmark, followed by the main event—the nuclear football analysis.

## Prerequisites

This book assumes familiarity with:

- **Calculus**: Integration, differential equations
- **Probability**: Random variables, probability distributions, expected values
- **Linear algebra**: Vectors, matrices, basic operations
- **Basic physics**: Energy, momentum, classical mechanics

No prior knowledge of nuclear physics is required—we'll build everything from first principles.

## A Note on Units

Nuclear physics has its own unit conventions. Throughout this book:

- **Energy**: Electron-volts (eV), with 1 MeV = 10⁶ eV common for fast neutrons
- **Cross sections**: Barns (b), where 1 barn = 10⁻²⁴ cm²
- **Length**: Centimeters (cm) for macroscopic distances, femtometers (fm = 10⁻¹³ cm) for nuclear scales
- **Mass**: Grams (g) or kilograms (kg) for bulk material, atomic mass units (u) for nuclides

The barn is a famously whimsical unit—nuclear physicists working on the Manhattan Project chose it because nuclear cross sections are "as big as a barn" compared to geometric nuclear sizes (which are ~10⁻²⁴ cm²).

## Let's Begin

Turn the page to start our journey into the atomic nucleus.
