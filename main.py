import math
from pprint import pprint
import time

import numpy as np
import matplotlib.pyplot as plt

m0 = 4 * math.pi * 10 ** -7

def get_coords(x0=-13.2 * 10 ** -2, y0=0, r=4 * 10 ** -2, m=0.63, count_of_magnets=3):
    delta_angle = 2 * math.pi / count_of_magnets
    list_of_main_coords = []
    list_of_second_coords = []
    magnet_moments = []
    for i in range(360):
        curr_angle = i * math.pi / 180
        coords_x_main = []
        coords_y_main = []
        coords_x_second = []
        coords_y_second = []

        magnet_vector_x = []
        magnet_vector_y = []

        for j in range(count_of_magnets):
            coords_x_main.append(r * math.cos(curr_angle + j * delta_angle))
            coords_y_main.append(r * math.sin(curr_angle + j * delta_angle))
            coords_x_second.append(x0 + r * math.cos(curr_angle + j * delta_angle))
            coords_y_second.append(y0 + r * math.sin(curr_angle + j * delta_angle))

            magnet_vector_x.append(m * math.cos(curr_angle + j * delta_angle))
            magnet_vector_y.append(m * math.sin(curr_angle + j * delta_angle))
        list_of_main_coords.append([coords_x_main, coords_y_main])
        list_of_second_coords.append([coords_x_second, coords_y_second])
        magnet_moments.append([magnet_vector_x, magnet_vector_y])
    return  list_of_main_coords, list_of_second_coords



def energy_matrix(x0=-13.2 * 10 ** -2, y0=0, r=4 * 10 ** -2, m=0.63, count_of_magnets=3):
    delta_angle = 2 * math.pi / count_of_magnets
    list_of_main_coords = []
    list_of_second_coords = []
    magnet_moments = []
    for i in range(360):
        curr_angle = i * math.pi / 180
        coords_x_main = []
        coords_y_main = []
        coords_x_second = []
        coords_y_second = []

        magnet_vector_x = []
        magnet_vector_y = []

        for j in range(count_of_magnets):
            coords_x_main.append(r * math.cos(curr_angle + j * delta_angle))
            coords_y_main.append(r * math.sin(curr_angle + j * delta_angle))
            coords_x_second.append(x0 + r * math.cos(curr_angle + j * delta_angle))
            coords_y_second.append(y0 + r * math.sin(curr_angle + j * delta_angle))

            magnet_vector_x.append(m * math.cos(curr_angle + j * delta_angle))
            magnet_vector_y.append(m * math.sin(curr_angle + j * delta_angle))
        list_of_main_coords.append([coords_x_main, coords_y_main])
        list_of_second_coords.append([coords_x_second, coords_y_second])
        magnet_moments.append([magnet_vector_x, magnet_vector_y])
    En_matrix = []
    for alpha in range(360):
        En_matrix.append([])
        for beta in range(360):
            Energy = 0
            for i in range(count_of_magnets):
                for j in range(count_of_magnets):
                    r = ((list_of_second_coords[beta][0][i] - list_of_main_coords[alpha][0][j]) ** 2 + (
                                list_of_second_coords[beta][1][i] - list_of_main_coords[alpha][1][j]) ** 2) ** 0.5
                    mr = ((list_of_second_coords[beta][0][i] - list_of_main_coords[alpha][0][j]) *
                          magnet_moments[alpha][0][j] + (
                                      list_of_second_coords[beta][1][i] - list_of_main_coords[alpha][1][j]) *
                          magnet_moments[alpha][1][j])
                    Bx = (m0 / (4 * math.pi)) * ((3 * mr) * (
                                list_of_second_coords[beta][0][i] - list_of_main_coords[alpha][0][j]) / r ** 5 -
                                                 magnet_moments[alpha][0][j] / r ** 3)
                    By = (m0 / (4 * math.pi)) * ((3 * mr) * (
                                list_of_second_coords[beta][1][i] - list_of_main_coords[alpha][1][j]) / r ** 5 -
                                                 magnet_moments[alpha][1][j] / r ** 3)
                    Energy += -(Bx * magnet_moments[beta][0][i] + By * magnet_moments[beta][1][i])
            En_matrix[alpha].append(Energy)
    return En_matrix


def interpolate_and_differentiate(matrix, a, b):
    # Преобразуем матрицу в numpy массив
    matrix = np.array(matrix)

    # Интерполируем значения по каждому из индексов
    interp_matrix = np.interp(np.linspace(0, 1, len(matrix[0])), np.arange(len(matrix[0])), matrix)

    # Вычисляем производную по второму индексу
    derivative_matrix = np.gradient(interp_matrix[a])[b]

    return interp_matrix




# e = energy_matrix(x0=11.5*10**-2)
# r = e
# # for i in range(3):
# #     for j in range(len(e)):
# #         r.append(e[j] + e[j] + e[j])
#
#
#
# a = open('energy2.csv', 'w')
# for x in r:
#     a.write(';'.join(map(str, x)) + '\n')
# a.close()
# exit()
#








import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1, )

a = list(map(lambda x: list(map(float, x.split(',')[:])), open('12.csv').readlines()))
k=3
c = 0
coords = get_coords(x0=11.5*10**-2, count_of_magnets=k)
colors = ['red', 'blue']
def animate(i):
    global c
    xs = []
    ys = []
    xs1 = []
    ys1 = []
    rxs, rys = [], []
    for i in range(2):
        for j in range(k):
            x, y = coords[i][int(a[c][i])][0][j], coords[i][int(a[c][i])][1][j]
            if i == 0:
                xs.append(0)
                ys.append(0)
                xs.append(x)
                ys.append(y)
                xs.append(0)
                ys.append(0)
            if i == 1:
                xs1.append(11.5 * 10 ** -2)
                ys1.append(0)
                xs1.append(x)
                ys1.append(y)
                xs1.append(11.5 * 10 ** -2)
                ys1.append(0)


    rxs.append(-0.1)
    rys.append(-0.1)
    rxs.append(0.17)
    rys.append(-0.1)
    rxs.append(-0.1)
    rys.append(0.1)
    rxs.append(0.17)
    rys.append(0.1)
    if c % 10 == 0:
        print(c, len(a))

    if c < len(a):
        c += 10
    ax1.clear()
    ax1.plot(xs1, ys1, c='red', linewidth=10)
    ax1.plot(xs, ys, c='blue', linewidth=10)
    ax1.plot(rxs, rys,'o', c='white')


ani = animation.FuncAnimation(fig, animate, repeat=False, interval=1)
# plt.show()
FFwriter = animation.PillowWriter()
ani.save('with.gif', writer=FFwriter)
#ani.save('orbita.gif', writer='imagemagick', fps=60)

