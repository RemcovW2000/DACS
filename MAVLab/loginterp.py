import numpy as np


def interpolate_log_scale(start, end, fraction):
    """
    Interpolates a value in a log scale between start and end given a fractional position.

    :param start: The starting value of the log scale (e.g., 10).
    :param end: The ending value of the log scale (e.g., 100).
    :param fraction: The fractional position between start and end (e.g., 0.61).
    :return: The interpolated value.
    """
    log_start = np.log10(start)
    log_end = np.log10(end)

    # Linear interpolation in log scale
    log_value = log_start + fraction * (log_end - log_start)

    # Convert back to the original scale
    value = 10 ** log_value

    return value


# Example usage
start_value = 1
end_value = 10
fraction = 12/59

interpolated_value = interpolate_log_scale(start_value, end_value, fraction)
print(f"The interpolated value is: {interpolated_value}")
