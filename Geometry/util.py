import math


def angle_difference(a, b):
	"""in radians"""
	diff = b - a
	while diff > math.pi:
		diff -= math.pi * 2
	while diff < -math.pi:
		diff += math.pi * 2
	return diff


def interpolate_angles(a, b, t):
	"""in radians"""
	return a + angle_difference(a, b) * t


def radians_to_degrees(radians: float) -> float:
	if radians is None:
		return None
	return radians / math.pi * 180.0


def opposite_angle(radians: float) -> float:
	result = radians - math.pi
	while result > math.pi:
		result -= math.pi * 2
	while result < -math.pi:
		result += math.pi * 2
	return result