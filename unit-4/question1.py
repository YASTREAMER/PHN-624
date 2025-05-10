import numpy as np
import matplotlib.pyplot as plt

def l_dot_s(l, j):
    """Calculates l.s for spin-orbit coupling."""
    if j == l + 1/2:
        return l/2
    elif j == l - 1/2:
        return -(l + 1)/2
    else:
        return 0

def l_square(l):
    """Calculates l(l+1)."""
    return l * (l + 1)

def energy_harmonic(n, l, hbar_w):
    """Calculates the harmonic oscillator energy."""
    N = 2 * (n - 1) + l
    return (N + 3/2) * hbar_w

def energy_l_squared(n, l, hbar_w, D):
    """Calculates the energy with l^2 term."""
    return energy_harmonic(n, l, hbar_w) + D * l_square(l)

def energy_total(n, l, j, hbar_w, C, D):
    """Calculates the total energy with l.s and l^2 terms."""
    return energy_l_squared(n, l, hbar_w, D) + C * l_dot_s(l, j)

def generate_shell_diagram(hbar_w=1, C=-0.1, D=-0.0225, n_max=3, l_max=2):
    """Generates the nuclear shell model energy level diagram."""

    energy_levels_harmonic = []
    energy_levels_l_squared = []
    energy_levels_total = []
    labels = []
    l_values = []
    n_values = []

    for n in range(1, n_max + 1):
        for l in range(0, l_max + 1):
            if 2 * (n - 1) + l >= n_max * 2:
                continue
            if l != 0:
                j_values = [l - 1/2, l + 1/2]
            else:
                j_values = [1/2]
            for j in j_values:
                energy_levels_harmonic.append(energy_harmonic(n, l, hbar_w))
                energy_levels_l_squared.append(energy_l_squared(n, l, hbar_w, D))
                energy_levels_total.append(energy_total(n, l, j, hbar_w, C, D))
                labels.append(f"n={n}, l={l}, j={j}")
                l_values.append(l)
                n_values.append(n)

    # Plotting
    unique_energies_harmonic = sorted(list(set(energy_levels_harmonic)))
    unique_energies_l_squared = sorted(list(set(energy_levels_l_squared)))
    unique_energies_total = sorted(list(set(energy_levels_total)))

    plt.figure(figsize=(10, 10))

    # Harmonic Oscillator Levels
    harmonic_positions = {}
    for e in unique_energies_harmonic:
        harmonic_positions[e] = e
        plt.hlines(y=e, xmin=0.5, xmax=1.5, color='orange', linestyle='-', linewidth=1)

    # l^2 Splitting and connection lines
    l_squared_positions = {}
    for e in unique_energies_l_squared:
        l_squared_positions[e] = e

    for e_harmonic, e_l_squared in zip(energy_levels_harmonic, energy_levels_l_squared):
        plt.plot([1.5, 2.5], [harmonic_positions[e_harmonic], l_squared_positions[e_l_squared]], color='black', linestyle='-', linewidth=0.5)
    
    for e in unique_energies_l_squared:
        plt.hlines(y=e, xmin=2.5, xmax=3.5, color='blue', linewidth=1)

    # Total Splitting (l.s) and connection lines
    total_positions = {}
    for e in unique_energies_total:
        total_positions[e] = e

    for e_l_squared, e_total in zip(energy_levels_l_squared, energy_levels_total):
        plt.plot([3.5, 4.5], [l_squared_positions[e_l_squared], total_positions[e_total]], color='black', linestyle='-', linewidth=0.5)

    for e in unique_energies_total:
        indices = [j for j, energy_val in enumerate(energy_levels_total) if energy_val == e]
        level_labels = [labels[j] for j in indices]
        label_text = ", ".join(level_labels)
        plt.text(4.6, total_positions[e], label_text, verticalalignment='center')

    # Updating colors for l.s levels
    for e in unique_energies_total:
        plt.hlines(y=e, xmin=4.5, xmax=5.5, color='red', linewidth=1)

    plt.yticks(list(total_positions.values()), [f"{e:.3f}" for e in unique_energies_total])
    plt.xticks([1, 3, 5], ["Harmonic", "l^2", "l.s"])
    plt.ylabel("Energy")
    plt.title("Nuclear Shell Model Energy Levels")
    plt.tight_layout()
    plt.show()

# Example usage
generate_shell_diagram()
