import math
import random
import time
import drawsvg as dw
import numpy as np

from moebius import Point, hyperbolic_distance, hyperbolic_geodesic, find_isometry_matrix


def octagon_error_function(vertex_coordinates):
    vertices = [Point(vertex_coordinates[2 * i], vertex_coordinates[2 * i + 1]) for i in range(8)]
    side_lengths = [hyperbolic_distance(vertices[i], vertices[(i + 1) % 8]) for i in range(8)]
    diagonal_lengths = [hyperbolic_distance(vertices[i], vertices[(i + 2) % 8]) for i in range(8)]

    result = 0
    indices = [0, 1, 4, 5]
    for i in indices:
        result += (side_lengths[i] - side_lengths[i + 2]) ** 2

    angle_sum = 0
    for i in range(8):
        a = side_lengths[i]
        b = side_lengths[(i + 1) % 8]
        c = diagonal_lengths[i]
        cos = (np.cosh(a) * np.cosh(b) - np.cosh(c)) / (np.sinh(a) * np.sinh(b))
        angle_sum += np.arccos(cos)
    result += (angle_sum - 2 * math.pi) ** 2
    return result


def approximate_gradient(coordinates):
    delta = 1e-06

    gradient = [0] * 16
    for i in range(16):
        coordinates_plus = [coordinates[j] if j != i else coordinates[j] + delta for j in range(16)]
        coordinates_minus = [coordinates[j] if j != i else coordinates[j] - delta for j in range(16)]
        gradient[i] = ((octagon_error_function(coordinates_plus) - octagon_error_function(coordinates_minus))
                       / (2 * delta))
    return np.array(gradient)


def descent(coordinates):
    step = 1e-03
    iterations = 200

    for i in range(iterations):
        h = step * approximate_gradient(coordinates)
        coordinates -= h
        if not check_coordinates(coordinates):
            return (False, coordinates)
    return (True, coordinates)


def check_coordinates(coordinates):
    magnitudes = np.array([np.hypot(coordinates[2 * i], coordinates[2 * i + 1]) for i in range(8)])
    return np.all(magnitudes < 1)


def draw_hyperbolic_segment(p1: Point, p2: Point, color):
    sx, sy = round(scale * p1.x) + offset, -round(scale * p1.y) + offset
    ex, ey = round(scale * p2.x) + offset, -round(scale * p2.y) + offset

    geodesic = hyperbolic_geodesic(p1, p2)

    if geodesic.is_line():
        p = dw.Path(stroke=color, fill='none', stroke_width=2)
        d.append(p.M(sx, sy).L(ex, ey))
    else:
        c = geodesic.center()
        theta1 = math.atan2(p1.y - c.y, p1.x - c.x)
        theta2 = math.atan2(p2.y - c.y, p2.x - c.x)
        r = round(geodesic.radius() * scale)

        from1to2 = theta2 - theta1 if theta2 - theta1 > 0 else theta2 - theta1 + 2 * math.pi
        from2to1 = theta1 - theta2 if theta1 - theta2 > 0 else theta1 - theta2 + 2 * math.pi
        sweep = (from1to2 > from2to1)

        p = dw.Path(stroke=color, fill='none', stroke_width=2)
        d.append(p.M(sx, sy).A(r, r, rot=0, large_arc=0, sweep=sweep, ex=ex, ey=ey))

start_time = time.time()

success_flag = False
octagon = []
attempt_counter = -1

while not success_flag:
    initial_position = []
    for j in range(8):
        angle = random.random() * (2 * math.pi / 8) + (2 * math.pi * j / 8)
        radius = random.random() * 0.3 + 0.65
        initial_position.extend([radius * math.cos(angle), radius * math.sin(angle)])
    success_flag, octagon = descent(initial_position)
    attempt_counter += 1

octagon_points = [Point(octagon[2 * i], octagon[2 * i + 1]) for i in range(8)]

scale = 375
offset = 400
d = dw.Drawing(offset * 2, offset * 2, id_prefix="pic")
d.append(dw.Circle(offset, offset, scale,
                   stroke='black', fill='none',
                   stroke_width=1))

colors = ["red", "green", "red", "green", "blue", "purple", "blue", "purple"]
for i in range(8):
    draw_hyperbolic_segment(octagon_points[i], octagon_points[(i + 1) % 8], colors[i])

paired_indices = [2, 3, 0, 1, 6, 7, 4, 5]
isometries = []
for j in range(8):
    A, B, C, D = (octagon_points[j],
                  octagon_points[(j + 1) % 8],
                  octagon_points[paired_indices[j]],
                  octagon_points[(paired_indices[j] + 1) % 8])
    isometries.append(find_isometry_matrix(A, B, D, C))

for j in range(8):
    isometry_matrix = np.eye(4)
    index = j
    for k in range(6):
        isometry_matrix = isometries[index] @ isometry_matrix
        index = (paired_indices[index] - 1) % 8
        for i in range(8):
            transformed_points = [octagon_points[i].transform(isometry_matrix) for i in range(8)]
            draw_hyperbolic_segment(transformed_points[i], transformed_points[(i + 1) % 8], colors[i])

d.save_svg('tiling.svg')

end_time = time.time()
elapsed = round(end_time - start_time, 2)
print(f"Successfully generated a tiling in {elapsed} seconds; failed attempts: {attempt_counter}.")
