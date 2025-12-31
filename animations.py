"""
Nuclear Physics Animation Suite - Implosion to Detonation
Uses actual simulation data from nuclear_physics.py and spatial_neutronics.py

Run: uv run manim -pql animations.py <SceneName>
     -p: preview  -ql: low quality (fast)  -qh: high quality
"""

from manim import *
import numpy as np

import nuclear_physics as nuc
import animation_data as anim_data


CORE_COLOR = RED
TAMPER_COLOR = ORANGE
PUSHER_COLOR = BLUE_C
HE_COLOR = GRAY
NEUTRON_COLOR = YELLOW


class ImplosionSequence(Scene):
    def construct(self):
        title = Text("Implosion Assembly", font_size=44).to_edge(UP)
        self.play(Write(title))

        data = anim_data.generate_implosion_data(n_points=200)

        scale = 2.5

        def radius_to_screen(r_m: float) -> float:
            return r_m * 100 * scale

        R_core_0 = radius_to_screen(data.R_core[0])
        R_tamper_0 = radius_to_screen(data.R_tamper[0])
        R_pusher_0 = radius_to_screen(data.R_pusher[0])
        R_outer_0 = min(radius_to_screen(data.R_outer[0]), 3.5)

        core = Circle(radius=R_core_0, color=CORE_COLOR, fill_opacity=0.8)
        tamper = Annulus(
            inner_radius=R_core_0,
            outer_radius=R_tamper_0,
            color=TAMPER_COLOR,
            fill_opacity=0.6,
        )
        pusher = Annulus(
            inner_radius=R_tamper_0,
            outer_radius=R_pusher_0,
            color=PUSHER_COLOR,
            fill_opacity=0.4,
        )
        he_shell = Annulus(
            inner_radius=R_pusher_0,
            outer_radius=R_outer_0,
            color=HE_COLOR,
            fill_opacity=0.3,
        )

        assembly = VGroup(he_shell, pusher, tamper, core)
        assembly.move_to(LEFT * 2.5)

        labels = (
            VGroup(
                Text("Pu-239", font_size=16, color=CORE_COLOR),
                Text("U-238", font_size=16, color=TAMPER_COLOR),
                Text("Al", font_size=16, color=PUSHER_COLOR),
                Text("HE", font_size=16, color=HE_COLOR),
            )
            .arrange(DOWN, aligned_edge=LEFT, buff=0.2)
            .to_edge(RIGHT)
            .shift(UP)
        )

        mass_labels = VGroup(
            Text(f"{data.core_mass:.1f} kg", font_size=14),
            Text(f"{data.tamper_mass:.0f} kg", font_size=14),
            Text(f"{data.pusher_mass:.0f} kg", font_size=14),
            Text(f"{data.outer_mass:.0f} kg", font_size=14),
        )
        for ml, l in zip(mass_labels, labels):
            ml.next_to(l, RIGHT, buff=0.3)

        self.play(Create(assembly), Write(labels), Write(mass_labels), run_time=2)

        comp_tracker = ValueTracker(0)

        comp_text = always_redraw(
            lambda: VGroup(
                MathTex(r"\rho/\rho_0 = "),
                DecimalNumber(
                    1.0 + (data.compression.max() - 1) * comp_tracker.get_value(),
                    num_decimal_places=2,
                ),
            )
            .arrange(RIGHT)
            .next_to(assembly, DOWN, buff=0.5)
        )

        time_text = always_redraw(
            lambda: VGroup(
                Text("t = ", font_size=20),
                DecimalNumber(
                    comp_tracker.get_value() * data.t[-1] * 1e6,
                    num_decimal_places=1,
                    unit=r"\mu s",
                ),
            )
            .arrange(RIGHT)
            .next_to(comp_text, DOWN)
        )

        self.add(comp_text, time_text)

        def update_assembly(mob):
            alpha = comp_tracker.get_value()
            idx = int(alpha * (len(data.t) - 1))

            R_core_new = radius_to_screen(data.R_core[idx])
            R_tamper_new = radius_to_screen(data.R_tamper[idx])
            R_pusher_new = radius_to_screen(data.R_pusher[idx])

            core.become(Circle(radius=R_core_new, color=CORE_COLOR, fill_opacity=0.8))
            tamper.become(
                Annulus(
                    inner_radius=R_core_new,
                    outer_radius=R_tamper_new,
                    color=TAMPER_COLOR,
                    fill_opacity=0.6,
                )
            )
            pusher.become(
                Annulus(
                    inner_radius=R_tamper_new,
                    outer_radius=R_pusher_new,
                    color=PUSHER_COLOR,
                    fill_opacity=0.4,
                )
            )
            he_shell.become(
                Annulus(
                    inner_radius=R_pusher_new,
                    outer_radius=R_outer_0,
                    color=HE_COLOR,
                    fill_opacity=0.3,
                )
            )
            VGroup(he_shell, pusher, tamper, core).move_to(LEFT * 2.5)

        assembly.add_updater(update_assembly)
        self.play(comp_tracker.animate.set_value(1.0), run_time=6, rate_func=smooth)
        assembly.remove_updater(update_assembly)

        crit_marker = Text("k = 1 CRITICAL", font_size=24, color=YELLOW)
        crit_marker.next_to(core, UP, buff=0.3)
        self.play(Write(crit_marker), core.animate.set_color(YELLOW))

        self.wait(1)

        supercrit = Text("SUPERCRITICAL", font_size=28, color=RED)
        supercrit.next_to(crit_marker, UP)
        self.play(
            Write(supercrit),
            core.animate.set_color(WHITE),
            Flash(core, color=WHITE, line_length=0.3),
        )

        self.wait(2)


class FVMGridVisualization(Scene):
    def construct(self):
        title = Text("Finite Volume Discretization", font_size=40).to_edge(UP)
        self.play(Write(title))

        subtitle = Text("2D Axisymmetric (r, z) Mesh", font_size=24, color=GRAY)
        subtitle.next_to(title, DOWN)
        self.play(Write(subtitle))

        fvm = anim_data.generate_fvm_cell_data(Nr=15, Nz=30)

        grid_width = 5.0
        grid_height = 4.0

        r_max = fvm.r_faces[-1]
        z_max = fvm.z_faces[-1]
        z_min = fvm.z_faces[0]

        def r_to_x(r):
            return (r / r_max) * grid_width - grid_width / 2

        def z_to_y(z):
            return ((z - z_min) / (z_max - z_min)) * grid_height - grid_height / 2

        grid_lines = VGroup()

        for r in fvm.r_faces[::2]:
            x = r_to_x(r)
            line = Line(
                start=np.array([x, z_to_y(z_min), 0]),
                end=np.array([x, z_to_y(z_max), 0]),
                stroke_width=0.5,
                color=GRAY,
            )
            grid_lines.add(line)

        for z in fvm.z_faces[::2]:
            y = z_to_y(z)
            line = Line(
                start=np.array([r_to_x(0), y, 0]),
                end=np.array([r_to_x(r_max), y, 0]),
                stroke_width=0.5,
                color=GRAY,
            )
            grid_lines.add(line)

        grid_lines.shift(LEFT * 2)
        self.play(Create(grid_lines), run_time=2)

        axis_r = Arrow(
            start=np.array([r_to_x(0) - 2, z_to_y(z_min) - 0.3, 0]),
            end=np.array([r_to_x(r_max) - 2 + 0.5, z_to_y(z_min) - 0.3, 0]),
            color=WHITE,
            buff=0,
        )
        axis_z = Arrow(
            start=np.array([r_to_x(0) - 2 - 0.3, z_to_y(z_min), 0]),
            end=np.array([r_to_x(0) - 2 - 0.3, z_to_y(z_max) + 0.5, 0]),
            color=WHITE,
            buff=0,
        )
        r_label = MathTex("r").next_to(axis_r, RIGHT)
        z_label = MathTex("z").next_to(axis_z, UP)

        self.play(Create(axis_r), Create(axis_z), Write(r_label), Write(z_label))

        phi_max = np.max(fvm.source)
        cells = VGroup()

        Nr, Nz = len(fvm.r_centers), len(fvm.z_centers)

        for i in range(0, Nr, 2):
            for j in range(0, Nz, 2):
                if fvm.D[i, j] > 0.01:
                    phi_norm = fvm.source[i, j] / phi_max if phi_max > 0 else 0
                    color = interpolate_color(BLUE, RED, min(phi_norm, 1.0))

                    r0, r1 = fvm.r_faces[i], fvm.r_faces[min(i + 2, Nr)]
                    z0, z1 = fvm.z_faces[j], fvm.z_faces[min(j + 2, Nz)]

                    x0, x1 = r_to_x(r0) - 2, r_to_x(r1) - 2
                    y0, y1 = z_to_y(z0), z_to_y(z1)

                    cell = Rectangle(
                        width=x1 - x0,
                        height=y1 - y0,
                        fill_color=color,
                        fill_opacity=0.7,
                        stroke_width=0,
                    ).move_to(np.array([(x0 + x1) / 2, (y0 + y1) / 2, 0]))
                    cells.add(cell)

        self.play(FadeIn(cells), run_time=2)

        equations = (
            VGroup(
                MathTex(
                    r"-\nabla \cdot (D \nabla \phi) + \Sigma_a \phi = \nu \Sigma_f \phi"
                ),
                MathTex(
                    r"\text{FVM: } \sum_{\text{faces}} D_f A_f \frac{\partial \phi}{\partial n} = V(\nu\Sigma_f - \Sigma_a)\phi"
                ),
            )
            .arrange(DOWN, buff=0.3)
            .scale(0.7)
            .to_edge(RIGHT)
            .shift(UP)
        )

        box = SurroundingRectangle(equations, color=BLUE, buff=0.2)
        self.play(Write(equations), Create(box))

        legend_title = Text("Fission Source", font_size=18)
        gradient = Rectangle(width=2, height=0.3, fill_opacity=1, stroke_width=0)
        gradient.set_color(color=BLUE)
        low_label = Text("Low", font_size=14).next_to(gradient, LEFT)
        high_label = Text("High", font_size=14).next_to(gradient, RIGHT)

        legend_title.next_to(gradient, UP)
        VGroup(legend_title, gradient, low_label, high_label).to_edge(DOWN)

        self.play(
            Write(legend_title), Create(gradient), Write(low_label), Write(high_label)
        )
        self.wait(3)


class FluxEvolutionScene(Scene):
    def construct(self):
        title = Text("Neutron Flux Evolution", font_size=40).to_edge(UP)
        self.play(Write(title))

        flux_data = anim_data.generate_flux_field(compression=2.5, Nr=40, Nz=80)

        width, height = 4.0, 6.0
        r_max = flux_data.R_max
        z_max = flux_data.Z_max / 2
        z_min = -z_max

        def create_flux_heatmap(phi: np.ndarray, sphere_r: float) -> VGroup:
            heatmap = VGroup()
            phi_max = np.max(phi)
            if phi_max <= 0:
                return heatmap

            Nr, Nz = phi.shape
            dr = r_max / Nr
            dz = (z_max - z_min) / Nz

            for i in range(0, Nr, 2):
                for j in range(0, Nz, 2):
                    r = flux_data.r[i]
                    z = flux_data.z[j]
                    dist = np.sqrt(r**2 + z**2)
                    if dist > sphere_r * 1.2:
                        continue

                    phi_val = phi[i, j]
                    if phi_val < 0.01 * phi_max:
                        continue

                    phi_norm = np.log10(phi_val / phi_max + 1e-10) / np.log10(1.0)
                    phi_norm = max(0, min(1, (phi_norm + 3) / 3))
                    color = interpolate_color(BLUE_E, RED, phi_norm)

                    x = (r / r_max) * width / 2
                    y = ((z - z_min) / (z_max - z_min) - 0.5) * height
                    cell_w = (dr / r_max) * width / 2 * 2
                    cell_h = (dz / (z_max - z_min)) * height * 2

                    rect = Rectangle(
                        width=cell_w,
                        height=cell_h,
                        fill_color=color,
                        fill_opacity=0.8,
                        stroke_width=0,
                    ).move_to(np.array([x, y, 0]))
                    heatmap.add(rect)

            return heatmap

        heatmap = create_flux_heatmap(flux_data.phi, flux_data.sphere_radius)
        heatmap.shift(LEFT * 2)

        sphere_screen_r = (flux_data.sphere_radius / r_max) * width / 2
        sphere_outline = Circle(radius=sphere_screen_r, color=WHITE, stroke_width=2)
        sphere_outline.shift(LEFT * 2)

        self.play(FadeIn(heatmap), Create(sphere_outline), run_time=2)

        info = (
            VGroup(
                MathTex(f"k_{{eff}} = {flux_data.k_eff:.4f}"),
                MathTex(f"\\rho/\\rho_0 = {flux_data.compression:.1f}"),
                MathTex(f"R = {flux_data.sphere_radius * 100:.2f}" + r"\text{ cm}"),
            )
            .arrange(DOWN, aligned_edge=LEFT)
            .scale(0.8)
            .to_edge(RIGHT)
            .shift(UP)
        )

        self.play(Write(info))

        profile_axes = (
            Axes(
                x_range=[0, 1.2, 0.2],
                y_range=[0, 1.2, 0.2],
                x_length=3,
                y_length=2.5,
                axis_config={"include_tip": False},
            )
            .to_edge(RIGHT)
            .shift(DOWN * 1.5)
        )

        profile_labels = profile_axes.get_axis_labels(
            x_label=MathTex(r"r/R"), y_label=MathTex(r"\phi/\phi_{max}")
        )

        self.play(Create(profile_axes), Write(profile_labels))

        mid_j = flux_data.Nz // 2
        r_norm = flux_data.r / flux_data.sphere_radius
        phi_profile = flux_data.phi[:, mid_j]
        phi_norm = (
            phi_profile / np.max(phi_profile)
            if np.max(phi_profile) > 0
            else phi_profile
        )

        valid = r_norm <= 1.2
        profile_curve = profile_axes.plot_line_graph(
            r_norm[valid],
            phi_norm[valid],
            add_vertex_dots=False,
            line_color=YELLOW,
            stroke_width=2,
        )

        self.play(Create(profile_curve))

        theory_label = MathTex(
            r"\phi(r) \propto \frac{\sin(\pi r/R')}{\pi r/R'}", font_size=28
        ).next_to(profile_axes, DOWN)
        self.play(Write(theory_label))

        self.wait(3)


class SupercriticalExcursion(Scene):
    def construct(self):
        title = Text("Supercritical Excursion", font_size=44).to_edge(UP)
        self.play(Write(title))

        exc = anim_data.generate_excursion_data()
        t_us = exc.t * 1e6

        sphere = Circle(radius=1.5, color=RED, fill_opacity=0.4)
        sphere.shift(LEFT * 3.5)
        self.play(Create(sphere))

        n_axes = Axes(
            x_range=[0, t_us[-1], t_us[-1] / 4],
            y_range=[0, 35, 5],
            x_length=4,
            y_length=2,
            axis_config={"include_tip": False},
        ).shift(RIGHT * 2.5 + UP * 2)
        n_label = Text("log₁₀(N)", font_size=16).next_to(n_axes, UP, buff=0.1)

        alpha_axes = Axes(
            x_range=[0, t_us[-1], t_us[-1] / 4],
            y_range=[-2, 10, 2],
            x_length=4,
            y_length=2,
            axis_config={"include_tip": False},
        ).shift(RIGHT * 2.5)
        alpha_label = Text("α (10⁸/s)", font_size=16).next_to(alpha_axes, UP, buff=0.1)

        T_axes = Axes(
            x_range=[0, t_us[-1], t_us[-1] / 4],
            y_range=[0, 12, 2],
            x_length=4,
            y_length=2,
            axis_config={"include_tip": False},
        ).shift(RIGHT * 2.5 + DOWN * 2)
        T_label = Text("T (10⁸ K)", font_size=16).next_to(T_axes, UP, buff=0.1)
        x_label = Text("t (μs)", font_size=14).next_to(T_axes, DOWN, buff=0.1)

        self.play(
            Create(n_axes),
            Write(n_label),
            Create(alpha_axes),
            Write(alpha_label),
            Create(T_axes),
            Write(T_label),
            Write(x_label),
        )

        time_tracker = ValueTracker(0)
        R_0 = exc.R[0]

        def get_current_index():
            t_current = time_tracker.get_value()
            idx = np.searchsorted(t_us, t_current)
            return min(idx, len(t_us) - 1)

        log_N = np.log10(exc.N + 1)
        log_N_scaled = log_N / log_N.max() * 30
        alpha_scaled = exc.alpha / 1e8
        T_scaled = exc.T / 1e8

        n_curve = always_redraw(
            lambda: n_axes.plot_line_graph(
                t_us[: get_current_index() + 1],
                log_N_scaled[: get_current_index() + 1],
                add_vertex_dots=False,
                line_color=YELLOW,
                stroke_width=2,
            )
            if get_current_index() > 0
            else VGroup()
        )

        alpha_curve = always_redraw(
            lambda: alpha_axes.plot_line_graph(
                t_us[: get_current_index() + 1],
                alpha_scaled[: get_current_index() + 1],
                add_vertex_dots=False,
                line_color=GREEN,
                stroke_width=2,
            )
            if get_current_index() > 0
            else VGroup()
        )

        T_curve = always_redraw(
            lambda: T_axes.plot_line_graph(
                t_us[: get_current_index() + 1],
                T_scaled[: get_current_index() + 1],
                add_vertex_dots=False,
                line_color=ORANGE,
                stroke_width=2,
            )
            if get_current_index() > 0
            else VGroup()
        )

        self.add(n_curve, alpha_curve, T_curve)

        def update_sphere(s):
            idx = get_current_index()
            R_ratio = exc.R[idx] / R_0
            new_radius = 1.5 * R_ratio
            T_frac = exc.T[idx] / exc.T.max()

            if T_frac < 0.3:
                color = interpolate_color(RED, ORANGE, T_frac / 0.3)
            elif T_frac < 0.7:
                color = interpolate_color(ORANGE, YELLOW, (T_frac - 0.3) / 0.4)
            else:
                color = interpolate_color(YELLOW, WHITE, (T_frac - 0.7) / 0.3)

            s.become(
                Circle(
                    radius=new_radius, color=color, fill_opacity=0.4 + 0.4 * T_frac
                ).move_to(LEFT * 3.5)
            )

        sphere.add_updater(update_sphere)

        yield_text = always_redraw(
            lambda: Text(
                f"Yield: {exc.yield_kt[get_current_index()]:.1f} kt", font_size=24
            ).next_to(sphere, DOWN, buff=0.3)
        )
        self.add(yield_text)

        self.play(
            time_tracker.animate.set_value(t_us[-1]), run_time=8, rate_func=linear
        )
        sphere.remove_updater(update_sphere)

        final_yield = Text(
            f"Final Yield: {exc.yield_kt[-1]:.1f} kt", font_size=32, color=YELLOW
        )
        final_yield.next_to(sphere, UP, buff=0.5)
        self.play(Write(final_yield), Flash(sphere, color=WHITE, line_length=0.5))

        self.wait(2)


class DisassemblyScene(Scene):
    def construct(self):
        title = Text("Hydrodynamic Disassembly", font_size=40).to_edge(UP)
        self.play(Write(title))

        exc = anim_data.generate_excursion_data()
        peak_idx = exc.i_peak_power
        start_idx = max(0, peak_idx - len(exc.t) // 4)

        t_slice = exc.t[start_idx:] - exc.t[start_idx]
        R_slice = exc.R[start_idx:]
        v_slice = exc.v[start_idx:]
        alpha_slice = exc.alpha[start_idx:]

        core = Circle(radius=1.5, color=WHITE, fill_opacity=0.6)
        core.shift(LEFT * 3)

        pressure_arrows = VGroup()
        n_arrows = 12
        for i in range(n_arrows):
            angle = i * 2 * PI / n_arrows
            arrow = Arrow(
                start=core.get_center(),
                end=core.get_center()
                + 0.8 * np.array([np.cos(angle), np.sin(angle), 0]),
                color=ORANGE,
                buff=0.3,
            )
            pressure_arrows.add(arrow)

        self.play(Create(core), Create(pressure_arrows))

        info_text = (
            VGroup(
                Text("Radiation Pressure:", font_size=20),
                MathTex(r"P = \frac{E_{rad}}{3V} \propto T^4"),
                Text("", font_size=12),
                Text("Acceleration:", font_size=20),
                MathTex(r"\frac{d^2R}{dt^2} = \frac{4\pi R^2 P}{M}"),
            )
            .arrange(DOWN, aligned_edge=LEFT)
            .to_edge(RIGHT)
            .shift(UP)
        )

        self.play(Write(info_text))

        v_axes = (
            Axes(
                x_range=[0, t_slice[-1] * 1e6, 0.2],
                y_range=[0, v_slice.max() / 1000 * 1.1, 100],
                x_length=4,
                y_length=2.5,
                axis_config={"include_tip": False},
            )
            .to_edge(RIGHT)
            .shift(DOWN * 1.5)
        )

        v_label = Text("Expansion v (km/s)", font_size=16).next_to(v_axes, UP)
        t_label = Text("t (μs)", font_size=14).next_to(v_axes, DOWN)

        self.play(Create(v_axes), Write(v_label), Write(t_label))

        time_tracker = ValueTracker(0)

        def get_idx():
            t_current = time_tracker.get_value() * 1e-6
            idx = np.searchsorted(t_slice, t_current)
            return min(idx, len(t_slice) - 1)

        v_curve = always_redraw(
            lambda: v_axes.plot_line_graph(
                t_slice[: get_idx() + 1] * 1e6,
                v_slice[: get_idx() + 1] / 1000,
                add_vertex_dots=False,
                line_color=BLUE,
                stroke_width=2,
            )
            if get_idx() > 0
            else VGroup()
        )

        self.add(v_curve)

        def update_core(c):
            idx = get_idx()
            scale = R_slice[idx] / R_slice[0]
            alpha_val = alpha_slice[idx]
            color = WHITE if alpha_val > 0 else BLUE
            c.become(
                Circle(
                    radius=min(1.5 * scale, 3.5),
                    color=color,
                    fill_opacity=max(0.1, 0.6 / scale),
                ).move_to(LEFT * 3)
            )

        def update_arrows(arrows):
            idx = get_idx()
            scale = R_slice[idx] / R_slice[0]
            for i, arrow in enumerate(arrows):
                angle = i * 2 * PI / n_arrows
                center = LEFT * 3
                r = 0.8 * scale
                arrow.become(
                    Arrow(
                        start=center,
                        end=center
                        + min(r, 2.5) * np.array([np.cos(angle), np.sin(angle), 0]),
                        color=interpolate_color(ORANGE, BLUE, min(1, (scale - 1) / 2)),
                        buff=0.3 * scale,
                    )
                )

        core.add_updater(update_core)
        pressure_arrows.add_updater(update_arrows)

        self.play(
            time_tracker.animate.set_value(t_slice[-1] * 1e6),
            run_time=6,
            rate_func=linear,
        )

        core.remove_updater(update_core)
        pressure_arrows.remove_updater(update_arrows)

        subcrit_text = Text("SUBCRITICAL (α < 0)", font_size=24, color=BLUE)
        subcrit_text.next_to(core, DOWN)
        self.play(Write(subcrit_text))

        self.wait(2)


class CompositeOverlay(Scene):
    def construct(self):
        title = Text("Nuclear Excursion - Multi-View", font_size=36).to_edge(UP)
        self.play(Write(title))

        exc = anim_data.generate_excursion_data()

        left_panel = Rectangle(width=5, height=5, color=WHITE, stroke_width=1)
        left_panel.shift(LEFT * 3.5 + DOWN * 0.5)
        right_panel = Rectangle(width=5, height=5, color=WHITE, stroke_width=1)
        right_panel.shift(RIGHT * 3.5 + DOWN * 0.5)

        self.play(Create(left_panel), Create(right_panel))

        left_title = Text("Core Cross-Section", font_size=18).next_to(left_panel, UP)
        right_title = Text("Time Evolution", font_size=18).next_to(right_panel, UP)
        self.play(Write(left_title), Write(right_title))

        core = Circle(radius=1.5, color=RED, fill_opacity=0.5)
        core.move_to(left_panel.get_center())

        inner_rings = VGroup()
        for r_frac in [0.25, 0.5, 0.75]:
            ring = Circle(
                radius=1.5 * r_frac, color=WHITE, stroke_width=0.5, stroke_opacity=0.3
            )
            ring.move_to(core.get_center())
            inner_rings.add(ring)

        self.play(Create(core), Create(inner_rings))

        t_us = exc.t * 1e6
        graph_center = right_panel.get_center()

        mini_axes = Axes(
            x_range=[0, t_us[-1], t_us[-1] / 2],
            y_range=[0, 1.2, 0.5],
            x_length=4,
            y_length=1.5,
            axis_config={"include_tip": False, "tick_size": 0.05},
        ).move_to(graph_center + UP * 1)

        mini_axes2 = Axes(
            x_range=[0, t_us[-1], t_us[-1] / 2],
            y_range=[-0.5, 1.2, 0.5],
            x_length=4,
            y_length=1.5,
            axis_config={"include_tip": False, "tick_size": 0.05},
        ).move_to(graph_center + DOWN * 1)

        label1 = Text("Yield (norm)", font_size=12).next_to(mini_axes, LEFT, buff=0.1)
        label2 = Text("α (norm)", font_size=12).next_to(mini_axes2, LEFT, buff=0.1)

        self.play(Create(mini_axes), Write(label1), Create(mini_axes2), Write(label2))

        yield_norm = exc.yield_kt / exc.yield_kt.max()
        alpha_norm = exc.alpha / exc.alpha.max()

        time_tracker = ValueTracker(0)

        def get_idx():
            t_val = time_tracker.get_value()
            return min(int(t_val / t_us[-1] * len(t_us)), len(t_us) - 1)

        yield_curve = always_redraw(
            lambda: mini_axes.plot_line_graph(
                t_us[: get_idx() + 1],
                yield_norm[: get_idx() + 1],
                add_vertex_dots=False,
                line_color=YELLOW,
                stroke_width=2,
            )
            if get_idx() > 0
            else VGroup()
        )

        alpha_curve = always_redraw(
            lambda: mini_axes2.plot_line_graph(
                t_us[: get_idx() + 1],
                alpha_norm[: get_idx() + 1],
                add_vertex_dots=False,
                line_color=GREEN,
                stroke_width=2,
            )
            if get_idx() > 0
            else VGroup()
        )

        self.add(yield_curve, alpha_curve)

        def update_core(c):
            idx = get_idx()
            R_scale = exc.R[idx] / exc.R[0]
            T_frac = exc.T[idx] / exc.T.max()
            color = interpolate_color(RED, WHITE, T_frac)
            new_r = min(1.5 * R_scale, 2.2)
            c.become(
                Circle(radius=new_r, color=color, fill_opacity=0.5).move_to(
                    left_panel.get_center()
                )
            )

        core.add_updater(update_core)

        time_display = always_redraw(
            lambda: Text(
                f"t = {time_tracker.get_value():.2f} μs", font_size=16
            ).to_edge(DOWN)
        )
        self.add(time_display)

        self.play(
            time_tracker.animate.set_value(t_us[-1]), run_time=8, rate_func=linear
        )
        core.remove_updater(update_core)

        self.wait(2)


class NarrativeSequence(Scene):
    def construct(self):
        title = Text("Nuclear Device Physics", font_size=44)
        subtitle = Text("From Implosion to Detonation", font_size=28, color=GRAY)
        VGroup(title, subtitle).arrange(DOWN).move_to(ORIGIN)

        self.play(Write(title), Write(subtitle))
        self.wait(1)
        self.play(FadeOut(title), FadeOut(subtitle))

        self.phase_implosion()
        self.phase_criticality()
        self.phase_excursion()
        self.phase_disassembly()
        self.final_summary()

    def phase_implosion(self):
        phase = Text("Phase 1: Implosion", font_size=36, color=BLUE).to_edge(UP)
        self.play(Write(phase))

        core = Circle(radius=2, color=CORE_COLOR, fill_opacity=0.6)
        shell = Annulus(
            inner_radius=2, outer_radius=2.8, color=HE_COLOR, fill_opacity=0.3
        )

        arrows = VGroup()
        for angle in np.linspace(0, 2 * PI, 8, endpoint=False):
            arr = Arrow(
                start=3.5 * np.array([np.cos(angle), np.sin(angle), 0]),
                end=2.5 * np.array([np.cos(angle), np.sin(angle), 0]),
                color=YELLOW,
                buff=0,
            )
            arrows.add(arr)

        self.play(Create(core), Create(shell), Create(arrows))
        self.play(core.animate.scale(0.6), shell.animate.scale(0.6), run_time=2)

        comp_text = MathTex(r"\rho \rightarrow 2.5\rho_0").next_to(core, DOWN)
        self.play(Write(comp_text))

        self.wait(1)
        self.play(FadeOut(VGroup(phase, core, shell, arrows, comp_text)))

    def phase_criticality(self):
        phase = Text("Phase 2: Criticality", font_size=36, color=YELLOW).to_edge(UP)
        self.play(Write(phase))

        core = Circle(radius=1.2, color=YELLOW, fill_opacity=0.6)
        self.play(Create(core))

        k_text = MathTex(r"k_{eff} = 1.0 \rightarrow \text{CRITICAL}").next_to(
            core, DOWN
        )
        self.play(Write(k_text))

        neutrons = VGroup()
        for _ in range(5):
            n = Dot(color=WHITE, radius=0.05)
            n.move_to(
                core.get_center() + 0.3 * np.random.randn(3) * np.array([1, 1, 0])
            )
            neutrons.add(n)

        self.play(FadeIn(neutrons))

        for _ in range(3):
            new_neutrons = VGroup()
            for n in neutrons[:10]:
                for _ in range(2):
                    new_n = Dot(color=WHITE, radius=0.05)
                    offset = 0.4 * np.random.randn(3) * np.array([1, 1, 0])
                    new_n.move_to(n.get_center() + offset)
                    if np.linalg.norm(new_n.get_center()) < 1.5:
                        new_neutrons.add(new_n)
            self.play(FadeOut(neutrons), FadeIn(new_neutrons), run_time=0.5)
            neutrons = new_neutrons

        supercrit = MathTex(r"k_{eff} > 1 \rightarrow \text{SUPERCRITICAL}", color=RED)
        supercrit.next_to(k_text, DOWN)
        self.play(Write(supercrit), core.animate.set_color(RED))

        self.wait(1)
        self.play(FadeOut(VGroup(phase, core, k_text, supercrit, neutrons)))

    def phase_excursion(self):
        phase = Text(
            "Phase 3: Supercritical Excursion", font_size=36, color=RED
        ).to_edge(UP)
        self.play(Write(phase))

        core = Circle(radius=1.2, color=RED, fill_opacity=0.6)
        self.play(Create(core))

        equations = (
            VGroup(
                MathTex(r"N(t) = N_0 e^{\alpha t}"),
                MathTex(r"\alpha \sim 10^8 \text{ s}^{-1}"),
                MathTex(r"\tau_{generation} \sim 10 \text{ ns}"),
            )
            .arrange(DOWN)
            .scale(0.8)
            .to_edge(RIGHT)
        )

        self.play(Write(equations))

        for i in range(5):
            color = interpolate_color(RED, WHITE, i / 5)
            self.play(core.animate.scale(1.1).set_color(color), run_time=0.3)

        energy = MathTex(r"E \sim 10^{13} \text{ J} \sim 15 \text{ kt TNT}")
        energy.next_to(core, DOWN)
        self.play(Write(energy), Flash(core, color=WHITE))

        self.wait(1)
        self.play(FadeOut(VGroup(phase, core, equations, energy)))

    def phase_disassembly(self):
        phase = Text("Phase 4: Disassembly", font_size=36, color=BLUE).to_edge(UP)
        self.play(Write(phase))

        core = Circle(radius=1.5, color=WHITE, fill_opacity=0.4)
        self.play(Create(core))

        text = (
            VGroup(
                Text("Radiation pressure:", font_size=20),
                MathTex(r"P \sim 10^{14} \text{ Pa (Gbar)}"),
                Text("Expansion terminates chain", font_size=20),
            )
            .arrange(DOWN)
            .to_edge(RIGHT)
        )

        self.play(Write(text))
        self.play(core.animate.scale(2.5).set_opacity(0.1), run_time=2)

        subcrit = Text("α < 0 : Subcritical", color=BLUE).next_to(core, DOWN)
        self.play(Write(subcrit))

        self.wait(1)
        self.play(FadeOut(VGroup(phase, core, text, subcrit)))

    def final_summary(self):
        summary = VGroup(
            Text("Summary", font_size=40),
            Text(""),
            Text("• Implosion: ~10 μs", font_size=24),
            Text("• Criticality → Supercritical", font_size=24),
            Text("• Excursion: ~1 μs", font_size=24),
            Text("• Peak: ~10¹³ J in ~100 ns", font_size=24),
            Text("• Disassembly terminates reaction", font_size=24),
        ).arrange(DOWN, aligned_edge=LEFT)

        self.play(Write(summary), run_time=3)
        self.wait(3)


if __name__ == "__main__":
    print("Nuclear Physics Animation Suite")
    print("=" * 50)
    print("\nRun with: uv run manim -pql animations.py <SceneName>")
    print("\nAvailable scenes:")
    print("  ImplosionSequence      - Mass shell compression")
    print("  FVMGridVisualization   - Finite volume mesh & flux")
    print("  FluxEvolutionScene     - 2D neutron flux heatmap")
    print("  SupercriticalExcursion - Full excursion with graphs")
    print("  DisassemblyScene       - Hydrodynamic expansion")
    print("  CompositeOverlay       - Multi-view simultaneous")
    print("  NarrativeSequence      - Complete physics story")
