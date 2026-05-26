import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# SECCIÓN 2: CAÍDA LIBRE (sin fricción)
# ============================================================

# Constantes
g = 9.8

# Condiciones iniciales
x0 = 0.0
v0 = 10.0

# Tiempo
tfinal = 10.0
dt = tfinal / 2000
time = np.arange(0, tfinal, dt)
lal = time.size

# Inicializar vectores numéricos
xx = np.zeros(lal)
vx = np.zeros(lal)
xx[0] = x0
vx[0] = v0

# Bucle de Euler (caída libre)
i = 1
while i < len(time) and xx[i-1] > -0.0001:
    ax = -g
    vx[i] = vx[i-1] + ax * dt
    xx[i] = xx[i-1] + vx[i-1] * dt
    i += 1

# Solución analítica (caída libre)
va = v0 - g * time          # v(t) = v0 - g·t
xa = x0 + v0*time - 0.5*g*time**2   # x(t) = x0 + v0·t - ½g·t²

# --- Plot 1: Velocidad numérica vs analítica (caída libre) ---
plt.figure(1)
plt.plot(time, vx, 'b', label='Velocidad numérica')
plt.plot(time, va, 'r--', label='Velocidad analítica')
plt.xlabel('Tiempo (s)', fontsize=14)
plt.ylabel('Velocidad (m/s)', fontsize=14)
plt.title('Sección 2 — Velocidad: numérica vs analítica (caída libre)')
plt.xlim(0, 2.1)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('plot1_velocidad_caida_libre.png', dpi=150)
plt.show()

# --- Plot 2: Posición numérica vs analítica (caída libre) ---
plt.figure(2)
plt.plot(time, xx, 'b', label='Posición numérica')
plt.plot(time, xa, 'r--', label='Posición analítica')
plt.xlabel('Tiempo (s)', fontsize=14)
plt.ylabel('Posición (m)', fontsize=14)
plt.title('Sección 2 — Posición: numérica vs analítica (caída libre)')
plt.xlim(0, 2.1)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('plot2_posicion_caida_libre.png', dpi=150)
plt.show()


# ============================================================
# SECCIÓN 3: CAÍDA CON FRICCIÓN Ff = -b·v
# ============================================================

# Parámetros
m = 1.0    # kg
b = 0.1    # kg/s
# Mismas condiciones iniciales: x0=0, v0=10 m/s

# Inicializar vectores numéricos con fricción
xx_f = np.zeros(lal)
vx_f = np.zeros(lal)
ax_f = np.zeros(lal)
xx_f[0] = x0
vx_f[0] = v0
ax_f[0] = (-m*g - b*v0) / m

# Bucle de Euler (con fricción)
i = 1
while i < len(time) and xx_f[i-1] > -0.0001:
    ax_f[i] = (-m*g - b*vx_f[i-1]) / m    # a = -g - (b/m)·v
    vx_f[i] = vx_f[i-1] + ax_f[i-1] * dt
    xx_f[i] = xx_f[i-1] + vx_f[i-1] * dt
    i += 1

# Solución analítica con fricción
# v(t) = (v0 + mg/b)·e^(-bt/m) - mg/b
# x(t) = -(m/b)·(v0 + mg/b)·e^(-bt/m) - (mg/b)·t + m/b·(v0 + mg/b)
# a(t) = -g - (b/m)·v(t)

va_f = (v0 + m*g/b) * np.exp(-b*time/m) - m*g/b
xa_f = -(m/b)*(v0 + m*g/b)*np.exp(-b*time/m) - (m*g/b)*time + (m/b)*(v0 + m*g/b)
aa_f = -g - (b/m)*va_f

# --- Plot 3: Aceleración numérica vs analítica (con fricción) ---
plt.figure(3)
plt.plot(time, ax_f, 'b', label='Aceleración numérica')
plt.plot(time, aa_f, 'r--', label='Aceleración analítica')
plt.xlabel('Tiempo (s)', fontsize=14)
plt.ylabel('Aceleración (m/s²)', fontsize=14)
plt.title('Sección 3 — Aceleración: numérica vs analítica (con fricción)')
plt.xlim(0, 2.5)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('plot3_aceleracion_friccion.png', dpi=150)
plt.show()

# --- Plot 4: Velocidad numérica vs analítica (con fricción) ---
plt.figure(4)
plt.plot(time, vx_f, 'b', label='Velocidad numérica')
plt.plot(time, va_f, 'r--', label='Velocidad analítica')
plt.xlabel('Tiempo (s)', fontsize=14)
plt.ylabel('Velocidad (m/s)', fontsize=14)
plt.title('Sección 3 — Velocidad: numérica vs analítica (con fricción)')
plt.xlim(0, 2.5)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('plot4_velocidad_friccion.png', dpi=150)
plt.show()

# --- Plot 5: Posición numérica vs analítica (con fricción) ---
plt.figure(5)
plt.plot(time, xx_f, 'b', label='Posición numérica')
plt.plot(time, xa_f, 'r--', label='Posición analítica')
plt.xlabel('Tiempo (s)', fontsize=14)
plt.ylabel('Posición (m)', fontsize=14)
plt.title('Sección 3 — Posición: numérica vs analítica (con fricción)')
plt.xlim(0, 2.5)
plt.legend(fontsize=12)
plt.grid(True)
plt.tight_layout()
plt.savefig('plot5_posicion_friccion.png', dpi=150)
plt.show()