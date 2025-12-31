"""
Manim Animations for Nuclear Physics Simulation
3Blue1Brown-style visualizations of neutronics and criticality

Run with: uv run manim -pql animations.py <SceneName>
  -p: preview after render
  -ql: low quality (fast), use -qm or -qh for higher quality
"""

from manim import *
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import jv  # Bessel functions

# Physics constants
BARN_TO_M2 = 1e-28
MEV_TO_J = 1.602176634e-13
AMU_TO_KG = 1.66053906660e-27
m_n = 1.674927471e-27  # neutron mass kg

# Pu-239 nuclear data (fission-spectrum averaged)
PU239_DATA = {
    "name": "Pu-239",
    "A": 239.0,
    "M": 239.0521634,
    "rho_0": 19.86e3,  # kg/m³
    "sigma_f": 1.800,  # barns
    "sigma_el": 4.566,
    "sigma_inel": 1.369,
    "sigma_c": 0.065,
    "nu_prompt": 3.165,
    "E_fission": 200.0,  # MeV
    "a_watt": 0.966,
    "b_watt": 2.842,
}


def get_neutron_velocity():
    """Average fission neutron velocity from Watt spectrum."""
    a, b = PU239_DATA["a_watt"], PU239_DATA["b_watt"]
    E_avg_MeV = 1.5 * a + 0.25 * a**2 * b
    E_avg_J = E_avg_MeV * MEV_TO_J
    return np.sqrt(2 * E_avg_J / m_n)


def compute_diffusion_params(compression=1.0):
    """Compute diffusion parameters for Pu-239."""
    data = PU239_DATA
    rho = data["rho_0"] * compression
    n = rho / (data["M"] * AMU_TO_KG)  # atoms/m³

    sigma_s = data["sigma_el"] + data["sigma_inel"]
    sigma_t = data["sigma_f"] + sigma_s + data["sigma_c"]
    sigma_a = data["sigma_f"] + data["sigma_c"]

    Sigma_f = data["sigma_f"] * BARN_TO_M2 * n
    Sigma_s = sigma_s * BARN_TO_M2 * n
    Sigma_a = sigma_a * BARN_TO_M2 * n
    Sigma_t = sigma_t * BARN_TO_M2 * n

    mu_0 = 2.0 / (3.0 * data["A"])
    Sigma_tr = Sigma_t - mu_0 * Sigma_s

    lambda_tr = 1.0 / Sigma_tr
    D = lambda_tr / 3.0

    c = Sigma_s / Sigma_t
    correction = (1.0 - c) / 5.0
    D_tr = D * (1.0 - correction)

    B_m_sq = (data["nu_prompt"] * Sigma_f - Sigma_a) / D_tr
    k_inf = data["nu_prompt"] * Sigma_f / Sigma_a

    v = get_neutron_velocity()
    tau = 1.0 / (v * Sigma_f * data["nu_prompt"])

    return {
        "D_tr": D_tr,
        "lambda_tr": lambda_tr,
        "B_m_sq": B_m_sq,
        "k_inf": k_inf,
        "tau": tau,
        "Sigma_f": Sigma_f,
        "Sigma_a": Sigma_a,
        "v": v,
        "c": c,
    }


def compute_critical_radius(compression=1.0):
    """Critical radius for bare Pu-239 sphere."""
    params = compute_diffusion_params(compression)
    c = params["c"]
    delta = 0.7104 * params["lambda_tr"] * (1.0 - 0.25 * (1.0 - c))
    R_prime_c = np.pi / np.sqrt(params["B_m_sq"])
    return R_prime_c - delta, params


def compute_critical_mass(compression=1.0):
    """Critical mass for bare Pu-239 sphere."""
    R_c, params = compute_critical_radius(compression)
    volume = (4.0 / 3.0) * np.pi * R_c**3
    mass = volume * PU239_DATA["rho_0"] * compression
    return mass, R_c, params


# =============================================================================
# Scene 1: Critical Mass Concept
# =============================================================================


class CriticalMassScene(Scene):
    """
    Visualize the concept of critical mass through neutron multiplication.
    Shows subcritical, critical, and supercritical regimes.
    """

    def construct(self):
        # Title
        title = Text("Critical Mass", font_size=48).to_edge(UP)
        self.play(Write(title))

        # Create sphere visualization
        sphere_radius = 1.5

        # Subcritical demonstration
        subcrit_label = Text("Subcritical (k < 1)", font_size=32, color=BLUE).next_to(
            title, DOWN
        )
        self.play(Write(subcrit_label))

        sphere = Circle(radius=sphere_radius, color=BLUE, fill_opacity=0.3)
        self.play(Create(sphere))

        # Animate neutrons dying out
        self.animate_neutron_chain(
            sphere, generations=4, k_eff=0.7, color=BLUE, start_neutrons=1
        )

        self.play(FadeOut(subcrit_label))

        # Critical demonstration
        crit_label = Text("Critical (k = 1)", font_size=32, color=YELLOW).next_to(
            title, DOWN
        )
        self.play(Write(crit_label), sphere.animate.set_color(YELLOW))

        self.animate_neutron_chain(
            sphere, generations=5, k_eff=1.0, color=YELLOW, start_neutrons=1
        )

        self.play(FadeOut(crit_label))

        # Supercritical demonstration
        supercrit_label = Text(
            "Supercritical (k > 1)", font_size=32, color=RED
        ).next_to(title, DOWN)
        self.play(Write(supercrit_label), sphere.animate.set_color(RED))

        self.animate_neutron_chain(
            sphere, generations=5, k_eff=1.8, color=RED, start_neutrons=1
        )

        # Show exponential growth equation
        equation = MathTex(
            r"N(t) = N_0 \cdot e^{\alpha t}", r"\quad \alpha = \frac{k-1}{\tau}"
        ).next_to(sphere, DOWN, buff=0.5)
        self.play(Write(equation))

        self.wait(2)

        # Clean up
        self.play(FadeOut(VGroup(title, supercrit_label, sphere, equation)))

        # Show critical mass calculation
        self.show_critical_mass_calculation()

    def animate_neutron_chain(
        self, sphere, generations, k_eff, color, start_neutrons=1
    ):
        """Animate neutron chain reaction."""
        neutrons = []
        center = sphere.get_center()
        radius = sphere.radius

        # Create initial neutron at center
        for _ in range(start_neutrons):
            n = Dot(center, color=color, radius=0.08)
            neutrons.append(n)
            self.play(FadeIn(n), run_time=0.3)

        for gen in range(generations):
            new_neutrons = []

            for neutron in neutrons:
                # Each neutron produces k_eff new neutrons on average
                # Use probabilistic rounding
                n_new = int(k_eff) + (1 if np.random.random() < (k_eff % 1) else 0)
                n_new = min(n_new, 4)  # Cap for visualization

                for _ in range(n_new):
                    # Random direction
                    theta = np.random.uniform(0, 2 * np.pi)
                    phi = np.random.uniform(0, np.pi)
                    r = np.random.uniform(0.3, 0.8) * radius

                    new_pos = center + np.array(
                        [
                            r * np.sin(phi) * np.cos(theta),
                            r * np.sin(phi) * np.sin(theta),
                            0,
                        ]
                    )

                    # Keep within sphere for visualization
                    new_n = Dot(new_pos, color=color, radius=0.08)
                    new_neutrons.append(new_n)

            # Animate: old neutrons fade, new ones appear
            if new_neutrons:
                self.play(
                    *[FadeOut(n) for n in neutrons],
                    *[FadeIn(n) for n in new_neutrons],
                    run_time=0.5,
                )
            else:
                self.play(*[FadeOut(n) for n in neutrons], run_time=0.5)

            neutrons = new_neutrons[:20]  # Limit for performance

            if not neutrons:
                break

        # Clean up remaining neutrons
        if neutrons:
            self.play(*[FadeOut(n) for n in neutrons], run_time=0.3)

    def show_critical_mass_calculation(self):
        """Show the physics of critical mass calculation."""
        title = Text("Pu-239 Critical Mass Calculation", font_size=40).to_edge(UP)
        self.play(Write(title))

        # Key equations
        equations = (
            VGroup(
                MathTex(
                    r"\text{Material buckling: } B_m^2 = \frac{\nu\Sigma_f - \Sigma_a}{D}"
                ),
                MathTex(
                    r"\text{Geometric buckling: } B_g^2 = \left(\frac{\pi}{R + \delta}\right)^2"
                ),
                MathTex(r"\text{Criticality: } B_m^2 = B_g^2"),
                MathTex(
                    r"\text{Critical radius: } R_c = \frac{\pi}{\sqrt{B_m^2}} - \delta"
                ),
            )
            .arrange(DOWN, aligned_edge=LEFT, buff=0.4)
            .next_to(title, DOWN, buff=0.5)
        )

        for eq in equations:
            self.play(Write(eq), run_time=1)

        self.wait(1)

        # Calculate and show actual values
        mass, R_c, params = compute_critical_mass(compression=1.0)

        results = (
            VGroup(
                Text(f"For bare Pu-239 sphere:", font_size=28),
                MathTex(f"R_c = {R_c * 100:.2f}" + r"\text{ cm}"),
                MathTex(f"M_c = {mass:.2f}" + r"\text{ kg}"),
                MathTex(f"k_\\infty = {params['k_inf']:.3f}"),
            )
            .arrange(DOWN, aligned_edge=LEFT, buff=0.3)
            .to_edge(RIGHT)
            .shift(DOWN)
        )

        box = SurroundingRectangle(results, color=YELLOW, buff=0.2)

        self.play(Write(results), Create(box))
        self.wait(3)

        self.play(FadeOut(VGroup(title, equations, results, box)))


# =============================================================================
# Scene 2: Neutron Flux Distribution
# =============================================================================


class NeutronFluxScene(Scene):
    """
    Visualize 2D neutron flux distribution in a sphere.
    Shows the fundamental mode (cosine/Bessel) shape.
    """

    def construct(self):
        title = Text("Neutron Flux Distribution", font_size=44).to_edge(UP)
        self.play(Write(title))

        # Compute critical parameters
        R_c, params = compute_critical_radius(compression=1.0)
        B = np.sqrt(params["B_m_sq"])

        # Create 2D flux visualization
        axes = Axes(
            x_range=[-1.5, 1.5, 0.5],
            y_range=[-1.5, 1.5, 0.5],
            x_length=5,
            y_length=5,
            axis_config={"include_tip": False},
        ).shift(LEFT * 2)

        axes_labels = axes.get_axis_labels(x_label="r/R_c", y_label="z/R_c")

        self.play(Create(axes), Write(axes_labels))

        # Create heatmap of flux
        # φ(r,z) ∝ j0(B*r) * cos(Bz) for cylinder, sin(Br)/r for sphere
        resolution = 50
        x_vals = np.linspace(-1.4, 1.4, resolution)
        y_vals = np.linspace(-1.4, 1.4, resolution)

        # Create flux field (spherical fundamental mode)
        flux_field = VGroup()
        max_flux = 1.0

        for i, x in enumerate(x_vals):
            for j, y in enumerate(y_vals):
                r_norm = np.sqrt(x**2 + y**2)
                if r_norm <= 1.0:
                    # Fundamental mode: sin(πr/R')/(r/R') normalized
                    if r_norm < 0.01:
                        flux = 1.0
                    else:
                        flux = np.sin(np.pi * r_norm) / (np.pi * r_norm)
                    flux = max(0, flux)
                else:
                    flux = 0

                if flux > 0.01:
                    # Map to color
                    color = interpolate_color(BLUE, RED, flux)
                    point = axes.c2p(x, y)
                    dot = Dot(point, radius=0.06, color=color, fill_opacity=flux)
                    flux_field.add(dot)

        self.play(FadeIn(flux_field), run_time=2)

        # Add sphere boundary
        sphere_boundary = Circle(radius=axes.x_length / 3, color=WHITE, stroke_width=2)
        sphere_boundary.move_to(axes.c2p(0, 0))
        self.play(Create(sphere_boundary))

        # Add 1D profile
        profile_axes = (
            Axes(
                x_range=[0, 1.2, 0.2],
                y_range=[0, 1.2, 0.2],
                x_length=4,
                y_length=3,
                axis_config={"include_tip": False},
            )
            .to_edge(RIGHT)
            .shift(UP * 0.5)
        )

        profile_labels = profile_axes.get_axis_labels(
            x_label=MathTex("r/R_c"), y_label=MathTex(r"\phi/\phi_0")
        )

        self.play(Create(profile_axes), Write(profile_labels))

        # Flux profile
        def flux_profile(r):
            if r < 0.01:
                return 1.0
            elif r <= 1.0:
                return np.sin(np.pi * r) / (np.pi * r)
            else:
                return 0

        flux_curve = profile_axes.plot(
            flux_profile, x_range=[0.01, 1.0], color=YELLOW, stroke_width=3
        )
        self.play(Create(flux_curve))

        # Add equation
        equation = MathTex(
            r"\phi(r) = \phi_0 \frac{\sin(\pi r / R')}{(\pi r / R')}"
        ).next_to(profile_axes, DOWN)
        self.play(Write(equation))

        # Highlight that flux → 0 at extrapolated boundary
        extrap_note = Text(
            "φ → 0 at extrapolated boundary R' = R + δ", font_size=20, color=GRAY
        ).next_to(equation, DOWN)
        self.play(Write(extrap_note))

        self.wait(3)
        self.play(
            FadeOut(
                VGroup(
                    title,
                    axes,
                    axes_labels,
                    flux_field,
                    sphere_boundary,
                    profile_axes,
                    profile_labels,
                    flux_curve,
                    equation,
                    extrap_note,
                )
            )
        )


# =============================================================================
# Scene 3: Supercritical Excursion
# =============================================================================


class SupercriticalExcursionScene(Scene):
    """
    Animate a supercritical excursion showing:
    1. Exponential neutron growth
    2. Energy deposition
    3. Hydrodynamic expansion (disassembly)
    4. Final yield
    """

    def construct(self):
        title = Text("Supercritical Excursion", font_size=44).to_edge(UP)
        self.play(Write(title))

        # Create visualization area
        sphere = Circle(radius=1.5, color=RED, fill_opacity=0.3, stroke_width=3)
        sphere_label = MathTex("R = 1.2 R_c").next_to(sphere, DOWN)
        self.play(Create(sphere), Write(sphere_label))

        # Create graphs for time evolution
        # Neutron population graph
        n_axes = (
            Axes(
                x_range=[0, 10, 2],
                y_range=[0, 30, 5],
                x_length=4,
                y_length=2.5,
                axis_config={"include_tip": False},
            )
            .to_edge(RIGHT)
            .shift(UP * 1.5)
        )
        n_labels = n_axes.get_axis_labels(
            x_label=MathTex(r"t/\tau"), y_label=MathTex(r"\log N")
        )

        # Temperature graph
        T_axes = (
            Axes(
                x_range=[0, 10, 2],
                y_range=[0, 10, 2],
                x_length=4,
                y_length=2.5,
                axis_config={"include_tip": False},
            )
            .to_edge(RIGHT)
            .shift(DOWN * 1.5)
        )
        T_labels = T_axes.get_axis_labels(
            x_label=MathTex(r"t/\tau"), y_label=MathTex(r"T (keV)")
        )

        self.play(Create(n_axes), Write(n_labels), Create(T_axes), Write(T_labels))

        # Simulate supercritical excursion
        # Simplified model: dN/dt = α*N, dE/dt = E_f * N * Σ_f
        R_ratio = 1.2  # R/R_c
        k_eff = R_ratio**2 * 1.0  # Approximate k scaling
        alpha = (k_eff - 1) / 1.0  # Rossi alpha (normalized)

        t_vals = np.linspace(0, 10, 100)
        N_vals = np.exp(alpha * t_vals)
        N_vals = np.minimum(N_vals, 1e12)  # Cap for log plot

        # Energy accumulation (integral of N)
        E_vals = np.cumsum(N_vals) * (t_vals[1] - t_vals[0])
        # Temperature ∝ E^(1/4) roughly (radiation dominated)
        T_vals = (E_vals / E_vals[-1]) ** 0.25 * 10  # Scale to keV

        # Animate the excursion
        n_tracker = ValueTracker(0)

        # Neutron curve (log scale representation)
        def get_n_curve():
            t_max = n_tracker.get_value()
            if t_max < 0.1:
                return VGroup()
            idx = int(t_max / 10 * len(t_vals))
            idx = max(1, min(idx, len(t_vals)))
            return n_axes.plot_line_graph(
                t_vals[:idx],
                np.log10(N_vals[:idx] + 1) * 3,  # Scale for visibility
                add_vertex_dots=False,
                line_color=YELLOW,
            )

        def get_T_curve():
            t_max = n_tracker.get_value()
            if t_max < 0.1:
                return VGroup()
            idx = int(t_max / 10 * len(t_vals))
            idx = max(1, min(idx, len(t_vals)))
            return T_axes.plot_line_graph(
                t_vals[:idx], T_vals[:idx], add_vertex_dots=False, line_color=ORANGE
            )

        n_curve = always_redraw(get_n_curve)
        T_curve = always_redraw(get_T_curve)

        self.add(n_curve, T_curve)

        # Sphere color and size change with temperature
        def sphere_updater(s):
            t = n_tracker.get_value()
            idx = int(t / 10 * len(T_vals))
            idx = max(0, min(idx, len(T_vals) - 1))
            temp_frac = T_vals[idx] / 10.0

            # Color: blue → red → white
            if temp_frac < 0.5:
                color = interpolate_color(RED, ORANGE, temp_frac * 2)
            else:
                color = interpolate_color(ORANGE, WHITE, (temp_frac - 0.5) * 2)

            # Size expansion (hydrodynamic disassembly)
            if t > 5:
                expansion = 1 + 0.5 * ((t - 5) / 5) ** 2
            else:
                expansion = 1.0

            s.set_fill(color, opacity=0.3 + 0.4 * temp_frac)
            s.set_stroke(color)
            s.scale_to_fit_width(3.0 * expansion)

        sphere.add_updater(sphere_updater)

        # Animate time evolution
        self.play(n_tracker.animate.set_value(10), run_time=5, rate_func=linear)
        sphere.remove_updater(sphere_updater)

        # Show final yield
        yield_text = Text("Yield ~ 15 kt", font_size=36, color=YELLOW)
        yield_text.next_to(sphere, UP)
        self.play(Write(yield_text))

        # Rossi alpha explanation
        rossi_eq = MathTex(
            r"\alpha = \frac{k_{eff} - 1}{\tau} \approx 10^8 \text{ s}^{-1}"
        ).to_edge(DOWN)
        self.play(Write(rossi_eq))

        self.wait(2)
        self.play(
            FadeOut(
                VGroup(
                    title,
                    sphere,
                    sphere_label,
                    n_axes,
                    n_labels,
                    T_axes,
                    T_labels,
                    n_curve,
                    T_curve,
                    yield_text,
                    rossi_eq,
                )
            )
        )


# =============================================================================
# Scene 4: Compression and Critical Mass
# =============================================================================


class CompressionScene(Scene):
    """
    Show how compression reduces critical mass.
    M_c ∝ 1/ρ² (for constant k_eff)
    """

    def construct(self):
        title = Text("Implosion: Compression Reduces Critical Mass", font_size=36)
        title.to_edge(UP)
        self.play(Write(title))

        # Create compression visualization
        # Initial sphere
        initial_radius = 2.0
        sphere = Circle(radius=initial_radius, color=BLUE, fill_opacity=0.3)

        # Compression arrows
        arrows = VGroup()
        n_arrows = 8
        for i in range(n_arrows):
            angle = i * 2 * np.pi / n_arrows
            start = np.array([2.5 * np.cos(angle), 2.5 * np.sin(angle), 0])
            end = np.array([1.8 * np.cos(angle), 1.8 * np.sin(angle), 0])
            arrow = Arrow(start, end, color=YELLOW, buff=0)
            arrows.add(arrow)

        self.play(Create(sphere), Create(arrows))

        # Create graph
        graph_axes = Axes(
            x_range=[1, 3, 0.5],
            y_range=[0, 10, 2],
            x_length=5,
            y_length=3,
            axis_config={"include_tip": False},
        ).to_edge(RIGHT)

        graph_labels = graph_axes.get_axis_labels(
            x_label=MathTex(r"\rho/\rho_0"), y_label=MathTex(r"M_c \text{ (kg)}")
        )

        self.play(Create(graph_axes), Write(graph_labels))

        # Plot M_c vs compression
        compressions = np.linspace(1.0, 3.0, 50)
        masses = []
        for c in compressions:
            m, _, _ = compute_critical_mass(compression=c)
            masses.append(m)

        mass_curve = graph_axes.plot_line_graph(
            compressions, masses, add_vertex_dots=False, line_color=GREEN
        )
        self.play(Create(mass_curve))

        # Animate compression
        compression_tracker = ValueTracker(1.0)

        def update_sphere(s):
            c = compression_tracker.get_value()
            # Radius scales as 1/c^(1/3)
            new_radius = initial_radius / (c ** (1 / 3))
            s.scale_to_fit_width(2 * new_radius)
            # Color intensity increases
            s.set_fill(opacity=0.3 * c)

        sphere.add_updater(update_sphere)

        # Compression indicator
        comp_text = always_redraw(
            lambda: MathTex(
                f"\\rho/\\rho_0 = {compression_tracker.get_value():.2f}"
            ).next_to(sphere, DOWN)
        )
        self.add(comp_text)

        # Dot on curve
        graph_dot = always_redraw(
            lambda: Dot(
                graph_axes.c2p(
                    compression_tracker.get_value(),
                    compute_critical_mass(compression_tracker.get_value())[0],
                ),
                color=RED,
            )
        )
        self.add(graph_dot)

        # Animate compression from 1x to 2.5x
        self.play(
            compression_tracker.animate.set_value(2.5), run_time=4, rate_func=smooth
        )

        sphere.remove_updater(update_sphere)

        # Show scaling law
        scaling_law = MathTex(r"M_c \propto \frac{1}{\rho^2}").to_edge(DOWN)
        self.play(Write(scaling_law))

        # Final values
        m_1, R_1, _ = compute_critical_mass(1.0)
        m_25, R_25, _ = compute_critical_mass(2.5)

        comparison = (
            VGroup(
                Text(f"At ρ₀: M_c = {m_1:.1f} kg", font_size=24),
                Text(f"At 2.5ρ₀: M_c = {m_25:.1f} kg", font_size=24),
            )
            .arrange(DOWN, aligned_edge=LEFT)
            .next_to(scaling_law, UP)
        )

        self.play(Write(comparison))
        self.wait(3)


# =============================================================================
# Scene 5: Fizzle Probability
# =============================================================================


class FizzleScene(Scene):
    """
    Visualize predetonation (fizzle) probability due to Pu-240.
    Shows Monte Carlo concept of spontaneous fission initiating chain.
    """

    def construct(self):
        title = Text("Predetonation Risk: The Fizzle Problem", font_size=36)
        title.to_edge(UP)
        self.play(Write(title))

        # Explanation
        explanation = (
            VGroup(
                Text("Pu-240 undergoes spontaneous fission", font_size=24),
                Text("→ Neutrons can start chain reaction early", font_size=24),
                Text("→ Reduced yield ('fizzle')", font_size=24),
            )
            .arrange(DOWN, aligned_edge=LEFT)
            .next_to(title, DOWN)
        )

        self.play(Write(explanation), run_time=2)
        self.wait(1)
        self.play(FadeOut(explanation))

        # Create compression timeline visualization
        timeline = Line(LEFT * 5, RIGHT * 5, color=WHITE)
        timeline.shift(UP * 1)

        # Mark key points
        start_dot = Dot(LEFT * 5 + UP, color=GREEN)
        crit_dot = Dot(LEFT * 2 + UP, color=YELLOW)
        peak_dot = Dot(RIGHT * 2 + UP, color=RED)

        start_label = Text("Start", font_size=20).next_to(start_dot, DOWN)
        crit_label = Text("k=1", font_size=20).next_to(crit_dot, DOWN)
        peak_label = Text("Peak ρ", font_size=20).next_to(peak_dot, DOWN)

        self.play(
            Create(timeline),
            Create(VGroup(start_dot, crit_dot, peak_dot)),
            Write(VGroup(start_label, crit_label, peak_label)),
        )

        # Supercritical window
        window = Rectangle(width=4, height=0.5, color=RED, fill_opacity=0.3).move_to(
            UP + RIGHT * 0
        )
        window_label = Text("Supercritical window", font_size=18, color=RED)
        window_label.next_to(window, UP)

        self.play(Create(window), Write(window_label))

        # Monte Carlo visualization
        mc_title = Text("Monte Carlo: 1000 trials", font_size=24).shift(DOWN * 0.5)
        self.play(Write(mc_title))

        # Simulate outcomes
        np.random.seed(42)
        n_trials = 1000
        pu240_fraction = 0.06  # 6% Pu-240

        # Spontaneous fission rate for Pu-240: ~1000 n/s/kg
        # At ~5 kg: ~5000 n/s
        # Supercritical window: ~10 μs
        # Probability of SF during window: ~0.05

        sf_rate = 1000 * pu240_fraction * 10  # n/s for our mass
        window_time = 10e-6  # 10 μs

        outcomes = []
        for _ in range(n_trials):
            # Poisson: number of SF events in window
            n_sf = np.random.poisson(sf_rate * window_time)
            if n_sf > 0:
                # Early initiation - check if it leads to fizzle
                # Probability of chain ≈ (k-1)/k during early supercritical
                p_chain = 0.3  # Average during ramp-up
                if np.random.random() < p_chain:
                    outcomes.append("fizzle")
                else:
                    outcomes.append("nominal")
            else:
                outcomes.append("nominal")

        n_fizzle = outcomes.count("fizzle")
        fizzle_rate = n_fizzle / n_trials * 100

        # Animated bar chart of results
        bar_nominal = Rectangle(
            width=2, height=3 * (1 - n_fizzle / n_trials), color=GREEN, fill_opacity=0.7
        ).shift(DOWN * 2 + LEFT * 2)

        bar_fizzle = Rectangle(
            width=2, height=3 * (n_fizzle / n_trials), color=RED, fill_opacity=0.7
        ).shift(DOWN * 2 + RIGHT * 2)

        bar_nominal.align_to(DOWN * 3.5, DOWN)
        bar_fizzle.align_to(DOWN * 3.5, DOWN)

        label_nominal = Text(f"Nominal\n{100 - fizzle_rate:.1f}%", font_size=20)
        label_nominal.next_to(bar_nominal, UP)

        label_fizzle = Text(f"Fizzle\n{fizzle_rate:.1f}%", font_size=20)
        label_fizzle.next_to(bar_fizzle, UP)

        self.play(
            GrowFromEdge(bar_nominal, DOWN),
            GrowFromEdge(bar_fizzle, DOWN),
            Write(label_nominal),
            Write(label_fizzle),
            run_time=2,
        )

        # Note about weapons-grade
        note = Text(
            f"At {pu240_fraction * 100:.0f}% Pu-240: ~{fizzle_rate:.0f}% fizzle risk",
            font_size=24,
            color=YELLOW,
        ).to_edge(DOWN)
        self.play(Write(note))

        self.wait(3)


# =============================================================================
# Scene 6: Geometry Comparison
# =============================================================================


class GeometryComparisonScene(Scene):
    """
    Compare critical masses for different geometries:
    - Sphere (optimal)
    - Cylinder
    - Ellipsoid
    """

    def construct(self):
        title = Text("Geometry Effects on Criticality", font_size=40)
        title.to_edge(UP)
        self.play(Write(title))

        subtitle = Text(
            "Surface-to-volume ratio determines neutron leakage",
            font_size=24,
            color=GRAY,
        ).next_to(title, DOWN)
        self.play(Write(subtitle))

        # Create three geometries
        sphere = Circle(radius=1, color=GREEN, fill_opacity=0.3)
        sphere_label = Text("Sphere", font_size=24).next_to(sphere, DOWN)

        # Cylinder (ellipse representation for 2D)
        cylinder = Ellipse(width=1.5, height=2.5, color=YELLOW, fill_opacity=0.3)
        cylinder_label = Text("Cylinder", font_size=24).next_to(cylinder, DOWN)

        # Ellipsoid (oblate)
        ellipsoid = Ellipse(width=2.5, height=1.2, color=RED, fill_opacity=0.3)
        ellipsoid_label = Text("Ellipsoid", font_size=24).next_to(ellipsoid, DOWN)

        shapes = (
            VGroup(
                VGroup(sphere, sphere_label),
                VGroup(cylinder, cylinder_label),
                VGroup(ellipsoid, ellipsoid_label),
            )
            .arrange(RIGHT, buff=1.5)
            .shift(UP * 0.5)
        )

        self.play(
            Create(sphere),
            Write(sphere_label),
            Create(cylinder),
            Write(cylinder_label),
            Create(ellipsoid),
            Write(ellipsoid_label),
        )

        # Surface to volume ratios
        # Sphere: S/V = 3/R (minimum for given V)
        # Cylinder (L=2R): S/V = 4/R
        # Oblate ellipsoid (a=2c): S/V ≈ 3.3/R_eq

        ratios = (
            VGroup(
                MathTex(r"\frac{S}{V} = \frac{3}{R}", color=GREEN),
                MathTex(r"\frac{S}{V} = \frac{4}{R}", color=YELLOW),
                MathTex(r"\frac{S}{V} \approx \frac{3.5}{R}", color=RED),
            )
            .arrange(RIGHT, buff=1.5)
            .shift(DOWN * 1.5)
        )

        self.play(Write(ratios))

        # Critical mass comparison
        masses = (
            VGroup(
                Text("M_c = 8.9 kg", font_size=24, color=GREEN),
                Text("M_c = 12.1 kg", font_size=24, color=YELLOW),
                Text("M_c = 10.5 kg", font_size=24, color=RED),
            )
            .arrange(RIGHT, buff=1.2)
            .shift(DOWN * 2.5)
        )

        self.play(Write(masses))

        # Key insight
        insight = Text(
            "Sphere minimizes surface area → minimum leakage → minimum critical mass",
            font_size=22,
            color=BLUE,
        ).to_edge(DOWN)
        self.play(Write(insight))

        self.wait(3)

        # Morph sphere to show optimality
        self.play(
            sphere.animate.scale(1.3),
            cylinder.animate.set_opacity(0.1),
            ellipsoid.animate.set_opacity(0.1),
        )

        optimal = Text("OPTIMAL", font_size=20, color=GREEN)
        optimal.next_to(sphere, UP)
        self.play(Write(optimal))

        self.wait(2)


# =============================================================================
# Scene 7: Complete Overview
# =============================================================================


class NuclearPhysicsOverview(Scene):
    """
    Complete overview combining all concepts.
    """

    def construct(self):
        # Title
        title = Text("Nuclear Criticality: Key Concepts", font_size=44)
        title.to_edge(UP)
        self.play(Write(title))

        # Create concept map
        concepts = VGroup(
            VGroup(
                Text("1. Neutron Multiplication", font_size=24, color=BLUE),
                MathTex(
                    r"k = \frac{\text{neutrons in gen } n+1}{\text{neutrons in gen } n}"
                ),
            ).arrange(DOWN, buff=0.2),
            VGroup(
                Text("2. Critical Mass", font_size=24, color=GREEN),
                MathTex(r"M_c = \frac{4\pi}{3} R_c^3 \rho"),
            ).arrange(DOWN, buff=0.2),
            VGroup(
                Text("3. Compression", font_size=24, color=YELLOW),
                MathTex(r"M_c \propto \rho^{-2}"),
            ).arrange(DOWN, buff=0.2),
            VGroup(
                Text("4. Rossi Alpha", font_size=24, color=ORANGE),
                MathTex(r"\alpha = \frac{k-1}{\tau} \sim 10^8 \text{ s}^{-1}"),
            ).arrange(DOWN, buff=0.2),
            VGroup(
                Text("5. Yield", font_size=24, color=RED),
                MathTex(r"Y = \epsilon \cdot M \cdot 17 \text{ kt/kg}"),
            ).arrange(DOWN, buff=0.2),
        )

        concepts.arrange_in_grid(rows=2, cols=3, buff=0.8)
        concepts.next_to(title, DOWN, buff=0.5)

        for concept in concepts:
            self.play(Write(concept), run_time=1)

        self.wait(2)

        # Final message
        final = Text(
            "Based on arXiv:1606.01670v1 (Aste 2016)", font_size=20, color=GRAY
        ).to_edge(DOWN)
        self.play(Write(final))

        self.wait(3)


# =============================================================================
# Main entry point
# =============================================================================

if __name__ == "__main__":
    print("Run with: uv run manim -pql animations.py <SceneName>")
    print("\nAvailable scenes:")
    print("  CriticalMassScene        - Neutron multiplication concept")
    print("  NeutronFluxScene         - 2D flux distribution")
    print("  SupercriticalExcursionScene - Time evolution of excursion")
    print("  CompressionScene         - Effect of compression on M_c")
    print("  FizzleScene              - Predetonation probability")
    print("  GeometryComparisonScene  - Shape effects on criticality")
    print("  NuclearPhysicsOverview   - Complete concept summary")
