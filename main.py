import numpy as np
import cmath
import h5py
from tqdm import tqdm

# 定义物理常数
EPSILON_0 = 8.854e-12  # 真空介电常数 (F/m)
MU_0 = 4 * np.pi * 1e-7  # 真空磁导率 (H/m)
C = 3e8  # 光速 (m/s)


def generate_random_parameters(n_samples: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate random structural and wave parameters."""
    # Parameter ranges
    x_range = (0.5, 2.0)  # mm
    t_range = (0.01, 0.1)  # mm
    d_range = (0.01, 0.05)  # mm
    h_range = (0.05, 0.2)  # mm

    x = np.random.uniform(*x_range, n_samples)
    t = np.random.uniform(*t_range, n_samples)
    d = np.random.uniform(*d_range, n_samples)
    h = np.random.uniform(*h_range, n_samples)

    return x, t, d, h


def calculate_epsilon_p(d, x, t, epsilon_a):
    """计算 epsilon_p"""
    ratio = d / (x - t)
    term1 = ratio * (2 - ratio) * epsilon_a
    term2 = (1 - ratio) ** 2
    return term1 + term2


def calculate_epsilon_r(d, x, t, epsilon_a):
    """计算 epsilon_r"""
    ratio = d / (x - t)
    term1 = (1 - ratio) ** 2 - ratio * (2 - ratio)
    term2 = (term1 * (1 - epsilon_a))
    term3 = cmath.sqrt(term1 ** 2 * (1 - epsilon_a) ** 2 + 4 * epsilon_a)
    return 0.5 * (term2 + term3)


def calculate_parallel(t, x, epsilon_b, epsilon_p):
    """计算 epsilon_parallel 或 mu_parallel"""
    ratio = t / x
    term1 = ratio * (2 - ratio) * epsilon_b
    term2 = (1 - ratio) ** 2 * epsilon_p
    return term1 + term2


def calculate_perpendicular(t, x, epsilon_r, epsilon_b):
    """计算 epsilon_perpendicular 或 mu_perpendicular"""
    ratio = t / x
    term1 = (1 - ratio) ** 2 - ratio * (2 - ratio)
    term2 = term1 * (epsilon_r - epsilon_b)
    term3 = cmath.sqrt(term1 ** 2 * (epsilon_r - epsilon_b) ** 2 + 4 * epsilon_r * epsilon_b)
    return 0.5 * (term2 + term3)


def calculate_k0(f):
    """
    计算波数 k0
    f: 频率 (Hz)
    返回: 波数 (rad/m)
    """
    omega = 2 * np.pi * f
    return omega * np.sqrt(EPSILON_0 * MU_0)


def calculate_reflection_coefficient(epsilon_perp, mu_perp, f, h):
    """
    计算反射系数
    epsilon_perp: 垂直方向介电常数
    mu_perp: 垂直方向磁导率
    f: 频率 (Hz)
    h: 厚度 (m)
    返回: 反射系数 (dB)
    """
    # 计算波数
    k0 = calculate_k0(f)

    # 计算复数平方根
    sqrt_epsilon = cmath.sqrt(epsilon_perp)
    sqrt_mu = cmath.sqrt(mu_perp)

    # 计算传播常数
    gamma = k0 * cmath.sqrt(mu_perp * epsilon_perp)

    # 计算指数项
    exp_term = cmath.exp(2j * gamma * h)

    # 计算分子和分母
    numerator = (sqrt_epsilon + sqrt_mu) + (sqrt_epsilon - sqrt_mu) * exp_term
    denominator = (sqrt_epsilon - sqrt_mu) + (sqrt_epsilon + sqrt_mu) * exp_term

    # 计算反射系数（dB）
    reflection = 20 * np.log10(abs(numerator / denominator))

    return reflection


def main():
    # 给定电磁参数
    epsilon_b = 4  # 骨架相对介电常数
    epsilon_a = complex(12, -1)  # 吸收层相对介电常数
    mu_b = 1  # 骨架相对磁导率
    mu_a = complex(3, -1)  # 吸收层相对磁导率

    # 给定结构参数
    # 蜂窝状结构的壁厚为2t，蜂窝孔之间的距离为2x
    # Generate frequency points
    f_hz = np.linspace(1e9, 10e9, 1000)
    combinations_per_freq = 1000000

    # Create HDF5 file and initialize datasets
    with h5py.File('honeycomb_reflection_freq_sweep.h5', 'w') as f:
        # Create groups for each type of data
        params = f.create_group('parameters')
        results = f.create_group('results')

        # Create datasets for structural parameters
        params.create_dataset('x', (len(f_hz), combinations_per_freq), dtype='float64')
        params.create_dataset('t', (len(f_hz), combinations_per_freq), dtype='float64')
        params.create_dataset('d', (len(f_hz), combinations_per_freq), dtype='float64')
        params.create_dataset('h', (len(f_hz), combinations_per_freq), dtype='float64')
        params.create_dataset('f_hz', data=f_hz)

        # Create datasets for results
        results.create_dataset('r', (len(f_hz), combinations_per_freq), dtype='float64')

        # Add metadata
        f.attrs['epsilon_b'] = epsilon_b
        f.attrs['epsilon_a_real'] = epsilon_a.real
        f.attrs['epsilon_a_imag'] = epsilon_a.imag
        f.attrs['mu_b'] = mu_b
        f.attrs['mu_a_real'] = mu_a.real
        f.attrs['mu_a_imag'] = mu_a.imag

        # Process each frequency point
        for i, freq in enumerate(tqdm(f_hz, desc="Processing frequencies")):
            # Generate random parameters for this frequency
            x, t, d, h = generate_random_parameters(combinations_per_freq)

            # Store parameters
            params['x'][i, :] = x
            params['t'][i, :] = t
            params['d'][i, :] = d
            params['h'][i, :] = h

            # Calculate reflection coefficients
            r = np.zeros(combinations_per_freq)

            for j in range(combinations_per_freq):
                epsilon_r = calculate_epsilon_r(d[j], x[j], t[j], epsilon_a)
                epsilon_perpendicular = calculate_perpendicular(t[j], x[j], epsilon_r, epsilon_b)

                mu_r = calculate_epsilon_r(d[j], x[j], t[j], mu_a)
                mu_perpendicular = calculate_perpendicular(t[j], x[j], mu_r, mu_b)

                r[j] = calculate_reflection_coefficient(epsilon_perpendicular, mu_perpendicular, freq, h[j])

            # Store results
            results['r'][i, :] = r

            # Flush data to disk periodically
            if i % 10 == 0:
                f.flush()


if __name__ == "__main__":
    main()