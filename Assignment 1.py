import numpy as np

A_years = np.array([[1928], [1932], [1936], [1948], [1952], [1956], [1960],
                    [1964], [1968], [1972], [1976], [1980], [1984], [1988],
                    [1992], [1996], [2000], [2004]])

A_years2 = np.square(A_years)


A_deg1 = np.array([[1, 1928], [1, 1932], [1, 1936], [1, 1948], [1, 1952], [1, 1956], [1, 1960],
                  [1, 1964], [1, 1968], [1, 1972], [1, 1976], [1, 1980], [1, 1984], [1, 1988],
                  [1, 1992], [1, 1996], [1, 2000], [1, 2004]])

A_deg2 = np.vstack([np.ones(len(A_years)), np.concatenate(A_years),
                    np.concatenate(np.power(A_years, 2))]).T

A_deg3 = np.vstack([np.ones(len(A_years)), np.concatenate(A_years),
                    np.concatenate(np.power(A_years, 2)), np.concatenate(np.power(A_years, 3))]).T

# print(a)

b = np.array([[12.2], [11.9], [11.5], [11.9], [11.5], [11.5], [11.0], [11.4], [11.0],
              [11.07], [11.08], [11.60], [10.97], [10.54], [10.82], [10.94], [10.75],
              [10.93]])

A_q1, A_r1 = np.linalg.qr(A_deg1)

print('A_q1 = \n', A_q1)
print('A_r1 = \n', A_r1)

A_q2, A_r2 = np.linalg.qr(A_deg2)

print('A_q2 = \n', A_q2)
print('A_r2 = \n', A_r2)

A_q3, A_r3 = np.linalg.qr(A_deg3)

print('A_q3 = \n', A_q3)
print('A_r3 = \n', A_r3)

lstsq_lin = np.linalg.lstsq(A_deg1, b, rcond=None)[0]

print('Numpy [x_1, x_2]T = \n', np.linalg.lstsq(A_deg1, b, rcond=None)[0])

print('Numpy [x_1, x_2, x_3]T = \n', np.linalg.lstsq(A_deg2, b, rcond=None)[0])
