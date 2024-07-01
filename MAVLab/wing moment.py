import numpy as np
import matplotlib.pyplot as plt


def elliptical_lift_distribution(y, L0, b):
    return L0 * np.sqrt(1 - (2 * y / b) ** 2)


def total_lift_to_L0(total_lift, b, n_points=1000):
    # Define a function to integrate the lift distribution for a given L0
    def integrate_lift(L0):
        y = np.linspace(0, b / 2, n_points)
        lift_distribution = elliptical_lift_distribution(y, L0, b)
        total_lift_calc = 2 * np.trapz(lift_distribution, y)
        return total_lift_calc

    # Use a numerical method to solve for L0
    from scipy.optimize import fsolve
    L0_initial_guess = 1
    L0_solution, = fsolve(lambda L0: integrate_lift(L0) - total_lift, L0_initial_guess)

    return L0_solution


def internal_moment(y, L0, b, n_points):
    # Compute the lift distribution over the half span
    y_half = np.linspace(0, b / 2, n_points)
    lift_distribution_half = elliptical_lift_distribution(y_half, L0, b)

    # Initialize the internal moment array
    moment_distribution_internal = np.zeros_like(y)

    # For each spanwise location, integrate the lift from y to half-span
    for i, yi in enumerate(y):
        # Lift distribution from yi to b/2
        y_segment = y_half[y_half >= yi]
        lift_segment = elliptical_lift_distribution(y_segment, L0, b)

        # Compute the moment arm (distance from yi to each segment point)
        moment_arm = y_segment - yi

        # Compute the internal moment
        moment_distribution_internal[i] = np.trapz(lift_segment * moment_arm, y_segment)

    return moment_distribution_internal
n_points = 1000

def shear_force(y, L0, b):
    # Shear force is the integral of the lift distribution from y to b/2
    shear_distribution = np.zeros_like(y)
    for i, yi in enumerate(y):
        y_segment = np.linspace(yi, b / 2, n_points - i)
        lift_segment = elliptical_lift_distribution(y_segment, L0, b)
        shear_distribution[i] = np.trapz(lift_segment, y_segment)
    return shear_distribution


def plot_lift_moment_shear(total_lift, b, n_points=1000):
    # Find L0 from the total lift
    L0 = total_lift_to_L0(total_lift, b, n_points)

    # Calculate distributions
    y = np.linspace(0, b / 2, n_points)
    lift_distribution = elliptical_lift_distribution(y, L0, b)
    moment_distribution_internal = internal_moment(y, L0, b, n_points)
    shear_distribution = shear_force(y, L0, b)

    # Plotting the results
    plt.figure(figsize=(12, 8))

    plt.subplot(3, 1, 1)
    plt.plot(y, lift_distribution, label='Lift Distribution')
    plt.xlabel('Spanwise location (y)')
    plt.ylabel('Lift per unit span (L(y))')
    plt.title('Elliptical Lift Distribution')
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(y, moment_distribution_internal, label='Internal Moment Distribution', color='orange')
    plt.axhline(0, color='black', linestyle='--', linewidth=1)
    plt.xlabel('Spanwise location (y)')
    plt.ylabel('Internal Moment (M(y))')
    plt.title('Internal Moment Curve through the Wing')
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.plot(y, shear_distribution, label='Shear Force Distribution', color='green')
    plt.axhline(0, color='black', linestyle='--', linewidth=1)
    plt.xlabel('Spanwise location (y)')
    plt.ylabel('Shear Force (V(y))')
    plt.title('Shear Force Distribution due to Lift')
    plt.legend()

    plt.tight_layout()
    plt.show()

    return L0

# Example usage
total_lift = 12.5*9.91*4  # Total lift in Newtons
b = 3  # wingspan in meters
L0_result = plot_lift_moment_shear(total_lift, b)
L0_result