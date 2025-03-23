import math
import numpy as np
import scipy as sp


output_decimals = 3
eps = 1e-06


class Point:
    x: float
    y: float

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return np.allclose([self.x, self.y], [other.x, other.y], atol=eps)

    def __repr__(self):
        x = round(self.x, output_decimals)
        y = round(self.y, output_decimals)
        return f"Point({{{x}, {y}}})"

    def create_circle(self, r=0.0):
        return GeneralizedCircle([1, -2 * self.x, -2 * self.y, self.x ** 2 + self.y ** 2 - r ** 2])

    def transform(self, isometry_matrix):
        return GeneralizedCircle(isometry_matrix @ self.create_circle().coefficients.T).center()


def dist_squared(p1: Point, p2: Point):
    return (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2


Q = np.array([[0, 0, 0, -1],
     [0, 1/2, 0, 0],
     [0, 0, 1/2, 0],
     [-1, 0, 0, 0]])


class GeneralizedCircle:
    def __init__(self, coefficients):
        if np.allclose(coefficients, [0.0, 0.0, 0.0, 0.0], atol=eps):
            raise ValueError("Equation coefficients too small in absolute value")
        if not math.isclose(coefficients[0], 0.0, abs_tol=eps):
            self.coefficients = np.array(coefficients, dtype=float) / coefficients[0]
        else:
            self.coefficients = np.array(coefficients, dtype=float)

    def __repr__(self):
        coefficients_rounded = np.round(self.coefficients, output_decimals)
        a, bx, by, c = tuple(coefficients_rounded)
        return f"circle {{{a}(x^2 + y^2) + {bx}x + {by}y + {c} == 0}}"

    def __eq__(self, other):
        return np.linalg.matrix_rank(np.stack((self.coefficients, other.coefficients))) == 1

    def __mul__(self, other):
        return self.coefficients @ Q @ other.coefficients.T

    def is_point(self):
        return math.isclose(self * self, 0.0, abs_tol=eps)

    def is_line(self):
        return math.isclose(self.coefficients[0], 0.0, abs_tol=eps)

    def center(self):
        if self.is_line():
            raise ValueError("The generalized circle is line and has no center")
        else:
            return Point(-self.coefficients[1] / (2 * self.coefficients[0]),
                         -self.coefficients[2] / (2 * self.coefficients[0]))

    def radius(self):
        if not self.is_line():
            return math.sqrt((self * self) / (2 * self.coefficients[0]))
        else:
            return np.inf

    def reflection_matrix(self):
        if self.is_point():
            raise ValueError("The reference circle has zero radius")
        return np.eye(4) - 2 * (np.outer(self.coefficients, self.coefficients) @ Q) / (self * self)

    def reflect_circle(self, circle):
        return GeneralizedCircle(self.reflection_matrix() @ circle.coefficients.T)

    def reflect_point(self, point: Point):
        return self.reflect_circle(point.create_circle()).center()


def orthogonal_circle(c1: GeneralizedCircle, c2: GeneralizedCircle, c3: GeneralizedCircle):
    equation_matrix = np.stack((c1.coefficients, c2.coefficients, c3.coefficients)) @ Q
    solution = sp.linalg.null_space(equation_matrix)
    return GeneralizedCircle(solution.T[0])


class CirclePencil:
    def __init__(self, c1: GeneralizedCircle, c2: GeneralizedCircle):
        self.coefficient_matrix = np.stack((c1.coefficients, c2.coefficients))
        self.circle1 = c1
        self.circle2 = c2

    def __repr__(self):
        return f"pencil through {self.circle1} and {self.circle2}"

    def orthogonal_pencil(self):
        equation_matrix = self.coefficient_matrix @ Q
        solution = sp.linalg.null_space(equation_matrix)
        d1 = GeneralizedCircle(solution.T[0])
        d2 = GeneralizedCircle(solution.T[1])
        return CirclePencil(d1, d2)


origin = Point(0, 0)
absolute = origin.create_circle(1)


def hyperbolic_geodesic(p1: Point, p2: Point):
    return orthogonal_circle(p1.create_circle(), p2.create_circle(), absolute)


def hyperbolic_midpoint_bisector(p1: Point, p2: Point):
    if p1 == p2:
        raise ValueError("Two given points coincide")
    pencil = CirclePencil(p1.create_circle(), p2.create_circle()).orthogonal_pencil()
    c1 = pencil.circle1
    c2 = pencil.circle2
    return orthogonal_circle(c1, c2, absolute)


def hyperbolic_distance(p1: Point, p2: Point):
    gamma = 2 * dist_squared(p1, p2) / ((1 - dist_squared(p1, origin)) * (1 - dist_squared(p2, origin)))
    return np.arccosh(1 + gamma)


def find_isometry_matrix(p1: Point, q1: Point, p2: Point, q2: Point):
    if not math.isclose(hyperbolic_distance(p1, q1), hyperbolic_distance(p2, q2), abs_tol=eps):
        raise ValueError("Two segments have different lengths")

    bisector1 = hyperbolic_midpoint_bisector(p1, p2)
    reflection_matrix1 = bisector1.reflection_matrix()
    q_prime = bisector1.reflect_point(q1)

    bisector2 = hyperbolic_midpoint_bisector(q_prime, q2)
    reflection_matrix2 = bisector2.reflection_matrix()

    return reflection_matrix2 @ reflection_matrix1
