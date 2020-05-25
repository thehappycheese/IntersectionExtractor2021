import math

class Vector2:
	"""Pure python Vector2 class"""
	
	def __init__(self, x: float = 0.0, y: float = 0.0):
		if type(x) is tuple:
			self.x = x[0]
			self.y = x[1]
		else:
			self.x = x
			self.y = y
	
	def copy(self):
		return Vector2(self.x, self.y)
	
	def magnitude_squared(self):
		return self.x*self.x + self.y*self.y
	
	def magnitude(self):
		return math.sqrt(self.x*self.x + self.y*self.y)
	
	def direction(self):
		"""
		:return: radians using math.atan2 (this is what ArcMap calls 'arithmetic' direction I think)
			if self.x == longitude and self.y == lattitude then 0° = East, 90° = North
		"""
		return math.atan2(self.y, self.x)
	
	def unit(self):
		ll = self.magnitude()
		if ll == 0:
			# what the user wants to happen here depends on the application. Here we simply hot-potato the problem by returning a zero vector
			# this allows us to be lazy by not checking for errors here or in the user code if we dont need to
			return Vector2(0.0, 0.0)
		return self / ll
	
	def scale_y(self, scale_y):
		result = self.copy()
		result.y *= scale_y
		return result
	
	def clip_y(self, minimum: float, maximum: float):
		result = self.copy()
		result.y = min(maximum, max(minimum, result.y))
		return result
	
	def round_y(self):
		result = self.copy()
		result.y = round(result.y/2)*2
		return result
	
	def __copy__(self):
		return Vector2(self.x, self.y)
	
	def __repr__(self):
		return f"Vector2({self.x:.2f}, {self.y:.2f})"
	
	def __str__(self):
		"""suitable for SVG"""
		return f"{self.x:.5f}".rstrip("0").rstrip(".") +f" {self.y:.5f}".rstrip("0").rstrip(".")
	
	def __add__(self, other):
		if type(other) is Vector2:
			return Vector2(self.x + other.x, self.y + other.y)
		return Vector2(self.x + other, self.y + other)
	
	def __sub__(self, other):
		if type(other) is Vector2:
			return Vector2(self.x - other.x, self.y - other.y)
		return Vector2(self.x - other, self.y - other)
	
	def __mul__(self, other):
		if type(other) == Vector2:
			raise TypeError("Unable to multiply Vector2 with Vector2. Please use .dot() or .cross() functions instead")
		return Vector2(self.x * other, self.y * other)
	
	def __truediv__(self, other):
		if type(other) == Vector2:
			raise TypeError("Unable to multiply Vector2 with Vector2. Please use .dot() or .cross() functions instead")
		return Vector2(self.x / other, self.y / other)
	
	def __neg__(self):
		return Vector2(-self.x, -self.y)
	
	def __pos__(self):
		return self.copy()