import math
from pprint import pprint
import time




m0 = 4 * math.pi * 10 ** -7

def get_coords(x0=-9 * 10 ** -2, y0=0, r=4 * 10 ** -2, m=0.63, count_of_magnets=3):
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



def energy_matrix(x0=-9.5 * 10 ** -2, y0=0, r=2.4 * 10 ** -2, m=0.63, count_of_magnets=3):
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

def energy_matrix_with_real(x0=-10 * 10 ** -2, y0=0, r=2.7 * 10 ** -2, m=0.63, count_of_magnets=3, x_resolution=6, y_resolution=6, x_size=0.7*10**-2, y_size=2.3*10**-2):

    m = m / (x_resolution * y_resolution)

    delta_angle = 2 * math.pi / count_of_magnets
    list_of_main_coords = []
    list_of_second_coords = []
    magnet_moments = []
    for i in range(120):
        curr_angle = i * math.pi / 180
        coords_x_main = []
        coords_y_main = []
        coords_x_second = []
        coords_y_second = []

        magnet_vector_x = []
        magnet_vector_y = []

        for j in range(count_of_magnets):
            coords_x_main.append([])
            coords_y_main.append([])
            coords_x_second.append([])
            coords_y_second.append([])
            magnet_vector_x.append([])
            magnet_vector_y.append([])
            for x_deps in range(x_resolution):
                coords_x_main[j].append([])
                coords_y_main[j].append([])
                coords_x_second[j].append([])
                coords_y_second[j].append([])
                magnet_vector_x[j].append([])
                magnet_vector_y[j].append([])

                for y_deps in range(y_resolution):
                    x_m = (x_deps + x_deps + 1 - x_resolution) * (x_size / (2 * x_resolution))
                    y_m = (y_deps + y_deps + 1 - y_resolution) * (y_size / (2 * y_resolution))
                    alpha = 0.000000001 + curr_angle + j * delta_angle
                    if int(alpha * 180 / math.pi) % 90 != 0 or (int(alpha * 180 / math.pi) // 90) % 2 == 0:
                        beta = math.atan(y_m / x_m)
                        x0_ = r * math.cos(curr_angle + j * delta_angle)
                        y0_ = r * math.sin(curr_angle + j * delta_angle)
                        if alpha != 0:
                            talpha = math.tan(alpha)
                        else:
                            talpha = math.tan(alpha)


                        tfi = math.tan(beta - math.pi / 2 + alpha)
                        if math.sin(alpha) != 0:
                            x1 = x_m / math.sin(alpha)
                        else:
                            x1 = 0

                        yc = talpha * x0_ - tfi * x0_

                        flag = (alpha * 180 / math.pi) // 90
                        # if flag == 1 or flag == 2 or flag == 5 or flag == 6 or flag == 9 or flag == 10:
                        #     x2 = (yc + talpha * x1) / (talpha - tfi)
                        #     y2 = talpha * (x2 - x1)
                        #     x2 = (yc + talpha * x1) / (talpha - tfi)
                        # else:
                        x2 = (yc + talpha * x1) / (talpha - tfi)
                        y2 = talpha * (x2 - x1)
                    else:
                        x2 = r * math.cos(curr_angle + j * delta_angle) + x_m
                        y2 = r * math.sin(curr_angle + j * delta_angle) + y_m

                    #print(x2, y2, alpha * 180 / math.pi, sep=';')
                    coords_x_main[j][x_deps].append(x2)
                    coords_y_main[j][x_deps].append(y2)
                    coords_x_second[j][x_deps].append(x0 + x2)
                    coords_y_second[j][x_deps].append(y0 + y2)

                    magnet_vector_x[j][x_deps].append(m * math.cos(curr_angle + j * delta_angle))
                    magnet_vector_y[j][x_deps].append(m * math.sin(curr_angle + j * delta_angle))
        list_of_main_coords.append([coords_x_main, coords_y_main])
        list_of_second_coords.append([coords_x_second, coords_y_second])
        magnet_moments.append([magnet_vector_x, magnet_vector_y])
    En_matrix = []
    for alpha in range(120):
        En_matrix.append([])
        for beta in range(120):
            Energy = 0
            for i in range(count_of_magnets):
                
                for j in range(count_of_magnets):
                    for x_deps_alpha in range(x_resolution):
                        for y_deps_alpha in range(y_resolution):
                            for x_deps_beta in range(x_resolution):
                                for y_deps_beta in range(y_resolution):
                                    r = ((list_of_second_coords[beta][0][i][x_deps_beta][y_deps_beta] - list_of_main_coords[alpha][0][j][x_deps_alpha][y_deps_alpha]) ** 2 + (
                                                list_of_second_coords[beta][1][i][x_deps_beta][y_deps_beta] - list_of_main_coords[alpha][1][j][x_deps_alpha][y_deps_alpha]) ** 2) ** 0.5
                                    mr = ((list_of_second_coords[beta][0][i][x_deps_beta][y_deps_beta] - list_of_main_coords[alpha][0][j][x_deps_alpha][y_deps_alpha]) *
                                          magnet_moments[alpha][0][j][x_deps_alpha][y_deps_alpha] + (
                                                      list_of_second_coords[beta][1][i][x_deps_beta][y_deps_beta] - list_of_main_coords[alpha][1][j][x_deps_alpha][y_deps_alpha]) *
                                          magnet_moments[alpha][1][j][x_deps_alpha][y_deps_alpha])
                                    Bx = (m0 / (4 * math.pi)) * ((3 * mr) * (
                                                list_of_second_coords[beta][0][i][x_deps_beta][y_deps_beta] - list_of_main_coords[alpha][0][j][x_deps_alpha][y_deps_alpha]) / r ** 5 -
                                                                 magnet_moments[alpha][0][j][x_deps_alpha][y_deps_alpha] / r ** 3)
                                    By = (m0 / (4 * math.pi)) * ((3 * mr) * (
                                                list_of_second_coords[beta][1][i][x_deps_beta][y_deps_beta] - list_of_main_coords[alpha][1][j][x_deps_alpha][y_deps_alpha]) / r ** 5 -
                                                                 magnet_moments[alpha][1][j][x_deps_alpha][y_deps_alpha] / r ** 3)
                                    Energy += -(Bx * magnet_moments[beta][0][i][x_deps_beta][y_deps_beta] + By * magnet_moments[beta][1][i][x_deps_beta][y_deps_beta])

            En_matrix[alpha].append(Energy)
        if alpha % 5 == 0:
            print(alpha)
    return En_matrix



def interpolate_and_differentiate(matrix, a, b):
    # Преобразуем матрицу в numpy массив
    matrix = np.array(matrix)

    # Интерполируем значения по каждому из индексов
    interp_matrix = np.interp(np.linspace(0, 1, len(matrix[0])), np.arange(len(matrix[0])), matrix)

    # Вычисляем производную по второму индексу
    derivative_matrix = np.gradient(interp_matrix[a])[b]

    return interp_matrix



start_time = time.time()
e = energy_matrix_with_real()
for i in range(len(e)):
    e[i] = e[i] * 3
e = e * 3
# for i in range(3):
#     for j in range(len(e)):
#         r.append(e[j] + e[j] + e[j])


a = open('8points_res.csv', 'w')
for x in e:
    a.write(';'.join(map(str, x)) + '\n')

print("--- %s seconds ---" % (time.time() - start_time))
exit()







import matplotlib.pyplot as plt
import matplotlib.animation as animation


fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1, )

a = list(map(lambda x: list(map(float, x.split(',')[:])), open('12.csv').readlines()))
k=4
c = 0
coords = get_coords(x0=-11*10**-2, count_of_magnets=k, r=4*10**-2)
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
                xs1.append(-13 * 10 ** -2)
                ys1.append(0)
                xs1.append(x)
                ys1.append(y)
                xs1.append(-13 * 10 ** -2)
                ys1.append(0)


    rxs.append(-0.2)
    rys.append(-0.1)
    rxs.append(0.1)
    rys.append(-0.1)
    rxs.append(-0.2)
    rys.append(0.1)
    rxs.append(0.1)
    rys.append(0.1)
    if c % 10 == 0:
        print(c, len(a))

    if c < len(a):
        c += 5
    ax1.clear()
    ax1.plot(xs1, ys1, c='red', linewidth=10)
    ax1.plot(xs, ys, c='blue', linewidth=10)
    ax1.plot(rxs, rys,'o', c='white')


ani = animation.FuncAnimation(fig, animate, repeat=False, interval=1)
# plt.show()
FFwriter = animation.PillowWriter()
ani.save('Калебки.gif', writer=FFwriter)

