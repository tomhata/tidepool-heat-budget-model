import json
import numpy as np
from numpy import pi, sin, cos, tan, arcsin, arccos, arctan
import pandas as pd


class TidePoolModel:
    def __init__(
        self,
        init_path="init.json",
        input_path="tenyear_final_serialdate.txt",
        auto_run=True,
        verbose=False,
    ):
        """Load default variable values from init.json file"""
        with open(init_path) as f:
            attr = json.load(f)
        for key, value in attr.items():
            setattr(self, key, value)
        self.verbose = verbose
        self.input_path = input_path
        self.import_data()
        self.remove_columns()
        self.fix_columns()
        self.add_columns()

        if auto_run:
            self.run_model()
    
    def set_time_interval(self):
        print('Current start time is:', self.start_time)
        self.start_time = input('New start time?')
        print('Current end time is:', self.end_time)
        self.end_time = input('New end time?')
    
    def run_model(self):
        self.df = self.df[self.start_time : self.end_time]
        self.set_init_values()
        self.calc_replenish()
        self.make_coeff_matrix()
        self.calc_solar_time()
        self.calc_irradiance()
        self.calc_transmittance()
        self.calc_shortwave()
        self.main_loop()

    def main_loop(self):
        """main calculation loop through time series"""
        for i in range(len(self.df)):
            # Set values for current time step
            self.idx = self.df.index[i]
            self.t_air = self.df["t_air"][i]
            self.wind_speed = self.df["wind_speed"][i]
            self.t_dist[-2] = self.df["t_ocean"][i]
            self.q_sw = self.df["q_sw"][i]

            # Calculate evaporative energy and mass flux if pool is not empty
            # and fill_pool is false. Use longwave emissivity of water.
            if (self.mass_water > 0) & (not self.fill_pool[i]):
                self.calc_evap()
                self.calc_longwave(i, self.emissivity_water)
            # If pool is dry or being refilled, evaporative flux is zero. If pool
            # is dry, use longwave emissivity of granite. If pool is not dry but
            # being refilled, longwave heat flux is set to zero anyway.
            else:
                self.q_evap = 0
                self.mass_flux = 0
                self.calc_longwave(i, self.emissivity_rock)
            # Refill pool with ocean temp seawater if fill_pool true.
            if self.fill_pool[i]:
                self.t_dist[0] = self.df["t_ocean"][i]
                self.q_lw = 0
                self.depth = self.depth_max
                self.r_surf = self.r_ini
                self.area_surf = pi * (self.r_surf) ** 2
                self.mass_water = self.mass_water_ini
                self.vol = self.vol_ini
                self.salinity = self.salinity_ini

            self.sub_loop(i)
            self.update_values()
            if self.verbose & (i % 1000 == 0):
                print("step", i, "out of", len(self.df))

    def sub_loop(self, i):
        """Subloop for finite element model using Runge-Kutta"""
        for j in range(self.sub_loop_steps):
            self.t_pool = self.t_dist[0]
            self.t_rock = self.t_dist[1]
            if (self.mass_water > 0) & (not self.fill_pool[i]):
                self.calc_cv_surf()
                self.calc_cv_bottom()
            elif self.fill_pool[i]:
                self.t_pool = self.df["t_ocean"][i]
                self.q_cv_surf = 0
                self.calc_cv_bottom()
            else:
                self.t_pool = self.t_rock
                self.t_dist[0] = self.t_rock
                self.calc_cv_surf()
                self.q_cv_bottom = 0
            self.update_coeff_matrix()
            self.rk4()

    def update_values(self):
        """Update parameters at the end of each loop of main loop."""
        self.t_pool = self.t_dist[0]
        self.t_rock = self.t_dist[1]
        self.mass_water = self.mass_water + self.mass_flux * self.time_interval
        if self.mass_water < (self.mass_water_min * self.mass_water_ini):
            self.mass_water = 0
        if self.mass_water > 0:
            self.salinity = self.mass_salt / self.mass_water * 1000
        if (self.salinity > self.salinity_max) or (self.mass_water == 0):
            self.salinity = self.salinity_max
        self.mass_pool = self.mass_water + self.mass_salt
        self.seawater_prop()
        self.update_pool()
        self.df.loc[self.idx, "q_lw"] = self.q_lw
        self.df.loc[self.idx, "q_cv_surf"] = self.q_cv_surf
        self.df.loc[self.idx, "q_cv_bottom"] = self.q_cv_bottom
        self.df.loc[self.idx, "q_evap"] = self.q_evap
        self.df.loc[self.idx, "r_surf"] = self.r_surf
        self.df.loc[self.idx, "t_pool"] = self.t_pool
        self.df.loc[self.idx, "depth"] = self.depth
        self.df.loc[self.idx, "volume"] = self.vol
        self.df.loc[self.idx, "mass_water"] = self.mass_water
        self.df.loc[self.idx, "mass_pool"] = self.mass_pool
        self.df.loc[self.idx, "salinity"] = self.salinity

    def import_data(self):
        """Import data from tab-delimited .txt file"""
        columns = [
            "year",
            "day_of_year",
            "hour:min",
            "irradiance_pyr",
            "t_air",
            "tide_height_pred",
            "tide_height",
            "wave_height",
            "wind_speed",
            "wind_dir",
            "t_ocean",
            "empty",
            "date",
        ]
        self.df = pd.read_csv(
            self.input_path, header=None, delimiter="\t", names=columns
        )

    def remove_columns(self):
        """Remove unnecessary columns"""
        columns_to_keep = [
            "day_of_year",
            "irradiance_pyr",
            "t_air",
            "tide_height",
            "wave_height",
            "wind_speed",
            "t_ocean",
        ]
        self.df = self.df[columns_to_keep]

    def fix_columns(self):
        """Adjust column values and add columns"""
        self.df["t_air"] += 273.15  # C to K
        self.df["t_ocean"] += 273.15  # C to K
        self.df["wave_height"] /= 100  # cm to m

    def add_columns(self):
        """Add columns and set index to datetime"""
        d_time = pd.to_datetime(
            range(0, len(self.df) * 10, 10),
            unit="m",
            origin=pd.Timestamp(self.start_time),
        )
        time_in_min = 60 * d_time.hour + d_time.minute
        self.df.insert(0, "d_time", d_time)
        self.df.insert(1, "time_in_min", time_in_min)
        self.df.set_index("d_time", inplace=True)
        cols_to_add = [
            "q_sw",
            "q_lw",
            "q_cv_surf",
            "q_cv_bottom",
            "q_evap",
            "r_surf",
            "t_pool",
            "depth",
            "volume",
            "mass_water",
            "mass_pool",
            "salinity",
        ]
        for col in cols_to_add:
            self.df[col] = np.nan

    def set_init_values(self):
        """Set initial values of model"""
        self.t_pool = self.df["t_ocean"][0]
        self.t_rock = self.t_pool
        self.t_dist = np.ones(self.num_nodes_rock + 2)
        self.t_dist[0:-1] = self.t_pool
        self.salinity = self.salinity_ini
        self.seawater_prop()

        self.area_bottom = pi * self.r_bottom ** 2
        self.sub_loop_steps = np.int(np.ceil(self.time_interval / self.dt))
        self.beta_const_rock = 1 / (
            self.density_rock * self.cp_rock * self.dist_node_rock * self.area_bottom
        )
        self.r_ini = self.r_bottom + self.depth_max / tan(self.angle_wall)
        self.vol_ini = (
            pi
            / 3
            * (self.r_bottom ** 2 + self.r_bottom * self.r_ini + self.r_ini ** 2)
            * self.depth_max
        )
        self.mass_pool_ini = self.density_water * self.vol_ini
        self.mass_salt = self.mass_pool_ini * self.salinity / (1000 + self.salinity)
        self.mass_water_ini = self.mass_pool_ini - self.mass_salt

        # Set current conditions to initial
        self.r_surf = self.r_ini
        self.area_surf = pi * self.r_surf ** 2
        self.mass_water = self.mass_water_ini
        self.mass_pool = self.mass_pool_ini
        self.vol = self.vol_ini
        self.depth = self.depth_max

    def seawater_prop(self):
        """
        Calculate and update seawater properties
        Inputs: t_pool (K), salinity (ppt)
        outputs:
            density_water : density(kg/m^3)  
            isobar_exp    : isobaric expansivity (K^-1) 
            latent_heat   : latent heat of evap (J/kg)
            conduct_water : thermal conductivity (W/(m*K))
            kinem_visc    : kinematic viscosity (m^2/s)
            cp_water      : specific heat (J/(kg*K))
            therm_diff    : thermal diffusivity (m^2/s)
        """

        T = self.t_pool - 273.15  # Temp in C
        # Set salinity to salinity_limit if a limit is set and is below current
        # salinity value.
        if (self.salinity_limit >= 0) & (self.salinity_limit < self.salinity):
            S = self.salinity_limit  # Salt in ppt
        else:
            S = self.salinity
        s = self.salinity / 1000  # Salt in fraction

        # |density_water| Density
        a = np.array(
            [
                9.9992293295e02,
                2.0341179217e-02,
                -6.1624591598e-03,
                2.2614664708e-05,
                -4.6570659168e-08,
            ]
        )
        b = np.array(
            [
                8.0200240891e02,
                -2.0005183488e00,
                1.6771024982e-02,
                -3.0600536746e-05,
                -1.6132224742e-05,
            ]
        )
        rho_w = a[0] + a[1] * T + a[2] * T ** 2 + a[3] * T ** 3 + a[4] * T ** 4
        D_rho = (
            b[0] * s
            + b[1] * s * T
            + b[2] * s * T ** 2
            + b[3] * s * T ** 3
            + b[4] * s ** 2 * T ** 2
        )
        self.density_water = rho_w + D_rho

        # |isobar_exp| Isobaric expansivity
        drho_wdT = a[1] + 2 * a[2] * T + 3 * a[3] * T ** 2 + 4 * a[4] * T ** 3
        dD_rhodT = (
            b[1] * s + 2 * b[2] * s * T + 3 * b[3] * s * T ** 2 + 2 * b[4] * s ** 2 * T
        )
        drho_sw_sharqdT = drho_wdT + dD_rhodT
        self.isobar_exp = -drho_sw_sharqdT / self.density_water

        # |latent_heat| latent heat of evaporation
        c = np.array(
            [
                2.5008991412e06,
                -2.3691806479e03,
                2.6776439436e-01,
                -8.1027544602e-03,
                -2.0799346624e-05,
            ]
        )
        hfg_w = c[0] + c[1] * T + c[2] * T ** 2 + c[3] * T ** 3 + c[4] * T ** 4
        self.latent_heat = hfg_w * (1 - 0.001 * S)

        # |k| thermal conductivity
        T_90 = 1.00024 * T  # convert from T_90 to T_68
        S_P = S / 1.00472  # convert from S to S_P
        self.conduct_water = 10 ** (
            np.log10(240 + 0.0002 * S_P)
            + 0.434
            * (2.3 - (343.5 + 0.037 * S_P) / (T_90 + 273.15))
            * (1 - (T_90 + 273.15) / (647.3 + 0.03 * S_P)) ** (1 / 3)
            - 3
        )

        # |kinem_visc| kinematic viscosity
        d = np.array(
            [
                1.5700386464e-01,
                6.4992620050e01,
                -9.1296496657e01,
                4.2844324477e-05,
                1.5409136040e00,
                1.9981117208e-02,
                -9.5203865864e-05,
                7.9739318223e00,
                -7.5614568881e-02,
                4.7237011074e-04,
            ]
        )
        mu_w = d[3] + 1 / (d[0] * (T + d[1]) ** 2 + d[2])
        viscA = d[4] + d[5] * T + d[6] * T ** 2
        viscB = d[7] + d[8] * T + d[9] * T ** 2

        # |mu| dynamic viscosity
        mu = mu_w * (1 + viscA * s + viscB * s ** 2)
        self.kinem_visc = mu / self.density_water

        # |cp_water| specific heat
        T68 = 1.00024 * (T + 273.15)  # convert from T_90 to T_68
        A = 5.328 - 9.76 * 10 ** (-2) * S + 4.04 * 10 ** (-4) * S ** 2
        B = -6.913 * 10 ** (-3) + 7.351 * 10 ** (-4) * S - 3.15 * 10 ** (-6) * S ** 2
        C = 9.6 * 10 ** (-6) - 1.927 * 10 ** (-6) * S + 8.23 * 10 ** (-9) * S ** 2
        D = 2.5 * 10 ** (-9) + 1.666 * 10 ** (-9) * S - 7.125 * 10 ** (-12) * S ** 2
        self.cp_water = 1000 * (A + B * T68 + C * (T68 ** 2) + D * (T68 ** 3))

        # |therm_diff| thermal diffusivity
        self.therm_diff = self.conduct_water / (self.density_water * self.cp_water)

    def calc_replenish(self):
        if self.replenish_on:
            amp_tide = (
                self.df["tide_height"] + self.swash_factor * self.df["wave_height"]
            )
            self.fill_pool = amp_tide > self.replenish_threshold

    def make_coeff_matrix(self):
        """
        Make matrix of coefficients to solve finite element model. Adds only
        conduction constants. Matrix must be updated with coefficients from
        other heat flux sources during each step through loop.
        """
        conduc_const = self.conductivity_rock / (
            self.density_rock * self.cp_rock * self.dist_node_rock ** 2
        )
        self.coeff_matrix = np.zeros((len(self.t_dist), len(self.t_dist)))
        # Conduction between surface rock node and adjacent node.
        self.coeff_matrix[1, 1] = -conduc_const
        self.coeff_matrix[1, 2] = conduc_const
        # Populate remaining matrix with conduction constants for bedrock.
        for u in range(2, (len(self.t_dist) - 2)):
            self.coeff_matrix[u, u - 1] = conduc_const
            self.coeff_matrix[u, u] = -2 * conduc_const
            self.coeff_matrix[u, u + 1] = conduc_const

    def calc_solar_time(self):
        """Reconcile solar and clock time. solar_time is a (n,) aray in minutes
        of the current day."""
        longitude = 121.88  # Longitude of HMS (degrees west)
        std_longitude = 120  # time zone longitude
        t_long = (self.std_longitude - self.longitude) * 4  # time correction (min)
        eot = 0.1237 * cos(
            ((2 * pi) / 365.242) * (self.df["day_of_year"] + 88.289)
        ) + 0.1654 * cos(((2 * pi) / 182.621) * (self.df["day_of_year"] - 127.029))
        self.solar_time = self.df["time_in_min"] + eot + t_long

    def calc_irradiance(self):
        """
        Calculate total solar irradiance from:
        zenith     : zenith angle of the sun
        latitude   : latitude (radians)
        omega      : hour angle
        delta      : declination angle
        solar_time : input, used to calculate omega
        """
        self.latitude = self.latitude * pi / 180  # deg to radians
        omega = 2 * pi * ((self.solar_time - (12 * 60)) / (24 * 60))
        delta = arcsin(0.4093 * cos(0.0172 * (self.df["day_of_year"] - 173)))
        self.zenith = arccos(
            sin(self.latitude) * sin(delta)
            + cos(self.latitude) * cos(delta) * cos(omega)
        )
        solar_altitude = (pi / 2) - self.zenith
        self.irradiance = self.df["irradiance_pyr"] / sin(solar_altitude)
        # Keep readings where solar altitude was above > 0.087 and solar wattage
        # was > 50
        self.irradiance[solar_altitude <= 0.087] = 0
        self.irradiance[self.df["irradiance_pyr"] <= 50] = 0
        # At extremely low angles, calculations give spuriously large irradiances
        self.irradiance[self.irradiance > 1200] = 1200
        self.irradiance[self.irradiance < 0] = 0
        # Assume daylight is between 5am and 7pm
        self.irradiance[self.df["time_in_min"] <= 300] = 0
        self.irradiance[self.df["time_in_min"] >= 1140] = 0
        self.zenith[self.df["time_in_min"] <= 300] = pi / 2
        self.zenith[self.df["time_in_min"] >= 1140] = pi / 2

    def calc_transmittance(self):
        """
        Calculate the fraction of light transmitted through water's surface
        to the rock substrate at bottom of pool. Angle of incidence is defined
        by (n,) array zenith. Reflectance of perpendicular and parallel components
        of polarized lights are calculated, and total reflectance is the average
        of the two values. Output variable transmit is the total transmittance of
        light [1-(total reflectance)] is a fraction between 0 and 1.
        """
        refrac_idx = 1.33  # refractive index of water
        reflect_perpend = (
            (sin(self.zenith) - tan(self.zenith) / refrac_idx)
            / (sin(self.zenith) + tan(self.zenith) / refrac_idx)
        ) ** 2
        reflect_parallel = (
            (tan(self.zenith) - sin(self.zenith) / refrac_idx)
            / (tan(self.zenith) + sin(self.zenith) / refrac_idx)
        ) ** 2
        reflect_tot = 0.5 * (reflect_perpend + reflect_parallel)
        self.transmit = 1 - reflect_tot

    def calc_shortwave(self):
        """
        Calculate energy input through shortwave radiation q_sw. q_sw is calculated as
        the product of solar irradiance (self.irradiance), transmitted through water
        (self.transmit), across a projected area of the pool bottom (area_proj), and
        absorbed by the rock substrated (self.absorptivity_rock). q_sw is directly
        inserted into dataframe, as it can be calculated for the entire time series at
        once.
        """
        area_proj = self.area_bottom * cos(self.zenith)
        q_sw = area_proj * self.irradiance * self.transmit * self.absorptivity_rock
        self.df["q_sw"] = q_sw

    def calc_longwave(self, i, emissivity_surf):
        """Calculate longwave emissivity"""
        # Air emissivity also depends on relative humidity, but this effect is not
        # significant. The equation is: 7 .37e-7 * t_air**2.46 * humid_rel**0.143.
        # Decreasing relative humidity from 1 to 0.75 only changes air emissivity
        # by 5%. Therefore it is currently omitted from calculation.
        if (self.df["time_in_min"][i] <= 300) or (self.df["time_in_min"][i] > 11400):
            emissivity_air = 1
        else:
            emissivity_air = 7.37e-7 * self.t_air ** 2.46
        # Calculate view factor if view_factor_on is true. Otherwise assume
        # view_factor to be 1 (full view of sky).
        if self.view_factor_on:
            view_factor = self.calc_view_factor()
        else:
            view_factor = 1

        self.q_lw = self.sb_const * emissivity_surf * self.area_surf * view_factor * self.t_air ** 4 * (
            emissivity_air - 1
        ) + 4 * self.sb_const * emissivity_surf * self.area_surf * view_factor * self.t_air ** 3 * (
            self.t_air - self.t_pool
        )

    def calc_view_factor(self):
        # Calculates the view factor from the surface of the pool, which is the
        # fraction of sky that is "viewable" by the surface. View factor is
        # used to calculate long wave emissivity. The calculation is taken from
        # the center of the surface of the pool.
        depth_from_top = self.depth_max - self.depth
        if depth_from_top == 0:
            view_factor = 1
        else:
            sky_angle = arctan(self.r_ini / depth_from_top)
            view_factor = sky_angle / (pi / 2)
        return view_factor

    def calc_evap(self):
        """Calculate evaporative energy and mass flux"""
        # t_film is the temperature of the air directly above water
        t_film = 0.5 * (self.t_air + self.t_pool)
        kinem_visc_air = -1.25e-5 + 9.28e-8 * t_film
        # binary diffusion coefficient between air and water.
        # regression data from Bolz and Tuve 1976. (m^2/s)
        binary_diff = -2.775e-6 + 4.479e-8 * t_film + 1.656e-10 * t_film ** 2
        # Reynolds number of air moving over pool. characteristic length = pool diam.
        reynolds = 2 * self.r_surf * self.wind_speed / kinem_visc_air
        # Schmidt number
        schmidt = kinem_visc_air / binary_diff
        # Sherwood number, convection over flat plate
        sherwood = 0.664 * reynolds ** (1 / 2) * schmidt ** (1 / 3)
        # h_diff, mass transfer coefficient (m/s)
        h_diff = sherwood * binary_diff / (2 * self.r_surf)
        sat_dens_air = self.calc_sat_density(self.t_air)
        sat_dens_pool = self.calc_sat_density(self.t_pool)
        humid_rel_surf = self.calc_humid_rel_surf()
        vapor_density_air = self.humid_rel * sat_dens_air
        vapor_density_surf = humid_rel_surf * sat_dens_pool
        # mass and energy flux calculation
        # mass_flux is the rate of evaporation of water (g/s)
        self.mass_flux = (
            h_diff * (vapor_density_air - vapor_density_surf) * self.area_surf / 1000
        )
        self.q_evap = self.latent_heat * self.mass_flux

    def calc_sat_density(self, temp):
        """
        4th order polynomial to calculate saturation density of water in air.
        """
        t_c = temp - 273.15
        coeffs = np.array([3.32e-6, 9.17e-5, 0.0121, 0.315, 4.92])
        sat_density = np.polyval(coeffs, t_c)
        return sat_density

    def calc_humid_rel_surf(self):
        """
        calculate relative humidity above a surface as a function of the
        surface's water content, as described in Luke Hunt's PhD thesis.
        """
        rel_water_content = 1 - self.salinity / 1000
        humid_rel_surf = 0.21 * np.log(rel_water_content) + 1.04
        if humid_rel_surf > 1:
            humid_rel_surf = 1
        elif humid_rel_surf < 0:
            humid_rel_surf = 0
        return humid_rel_surf

    def calc_cv_surf(self):
        """
        Calculates sensible heat convection across the surface of the
        tidepool due to forced convection of air flow across surface. 
        q_cv_surf is the rate of energy input or removed from the
        tidepool (W).
        """
        # Temperature of air film directly above water
        t_film = 0.5 * (self.t_air + self.t_pool)
        # Kinematic viscosity of air film (m^2 / s)
        kinem_visc_air = -1.25e-5 + 9.28e-8 * t_film
        # Thermal conductivity of air [W / (m*K)]
        conduct_air = 0.005013 + 7.2e-5 * t_film
        # Prandtl number of at at approx 300K.
        prandtl = 0.71
        # Reynolds number of air over tidepool. cl = diameter of pool.
        reynolds = 2 * self.r_surf * self.wind_speed / kinem_visc_air
        # Nusselt number of laminar flow over flat plate.
        nusselt = 0.664 * reynolds ** (1 / 2) * prandtl ** (1 / 3)
        # Convective heat transfer coefficient
        h_conv = nusselt * conduct_air / (2 * self.r_surf)
        self.q_cv_surf = h_conv * self.area_surf * (self.t_air - self.t_pool)

    def calc_cv_bottom(self):
        """
        Calculate heat transfer between tidepool and rock surface at the
        bottom of the pool via free (natural) convection. q_cv_bottom is
        the energy transfer rate (W) to the tidepool from the rock
        surface. Conversly, -(q_cv_bottom) is the energy transferred to
        the rock surface from the tidepool.
        """
        grav = 9.8
        # Characteristic length. For a flat plate, l_c = area / perimeter.
        # This value is subjective.
        l_c = self.r_bottom / 2
        # Rayleigh number. Positive if t_pool > t_rock.
        rayleigh = (
            grav
            * self.isobar_exp
            * np.absolute(self.t_pool - self.t_rock)
            * l_c ** 3
            / (self.therm_diff * self.kinem_visc)
        )
        # Calculation of nusselt number depends on whether t_pool or t_rock
        # is greater. See Incropera and DeWitt p.551.
        if self.t_pool == self.t_rock:
            nusselt = 0
        elif self.t_pool > self.t_rock:
            # Top surface of a cooled plate.
            nusselt = 0.27 * rayleigh ** (1 / 4)
        elif rayleigh < 1e7:
            # Top surface of a warmed plate. Equation depends on rayleigh num
            nusselt = 0.54 * rayleigh ** (1 / 4)
        else:
            nusselt = 0.15 * rayleigh ** (1 / 3)
        # heat transfer coefficient
        h_conv_free = self.conduct_water * nusselt / l_c
        # heat transfer from tidepool to rock.
        self.q_cv_bottom = -h_conv_free * self.area_bottom * (self.t_pool - self.t_rock)

    def update_coeff_matrix(self):
        """
        For use in matrix to solve rate of temperature change. Divide each heat
        flux source (excluding conduction) by this value to calculate total
        temperature change for Tidepool and layer of rock beneath. This needs to
        be recalculated each iteration, as water volume will change.
        """
        if self.mass_water > 0:
            beta_const_water = 1 / (
                self.cp_water * (self.mass_water * (1 + self.salinity / 1000))
            )
            # Add heat fluxes to tidepool node. Longwave radiation, evaporative
            # cooling, forced convection between air and tidepool, free convection
            # between tidepool and rock surface.
            if self.q_lw is None:
                print('q_lw is none')
            elif self.q_evap is None:
                print('q_evap is none')
            elif self.q_cv_surf is None:
                print('q_cv_surf is none')
            elif self.q_cv_bottom is None:
                print('q_cv_bottom is none')
            
            self.coeff_matrix[0, -1] = beta_const_water * (
                self.q_lw + self.q_evap + self.q_cv_surf + self.q_cv_bottom
            )
            # Add heat fluxes to surface rock node. Shortwave radiation and free
            # convection between tidepool and rock surface (-q_cv_bottom).
            self.coeff_matrix[1, -1] = self.beta_const_rock * (
                self.q_sw - self.q_cv_bottom
            )
        else:
            # If current water mass is 0, pool is dry. longwave radiation
            # and forced convection should therefore be removed from rock
            # surface.
            self.coeff_matrix[0, -1] = 0
            # Add heat fluxes to surface rock node. Shortwave radiation,
            # longwave radiation, and forced convection.
            self.coeff_matrix[1, -1] = self.beta_const_rock * (
                self.q_lw + self.q_sw + self.q_cv_surf
            )

    def rk4(self):
        """
        Resolve heat budget using 4th order Runge-Kutta
        y' = F(x,y)
        y0 = f(x0) -> yn = f(xn), initial values known
        time step h
        k1 = F( xn       , yn        )
        k2 = F( xn + h/2 , yn + k1/2 )
        k3 = F( xn + h/2 , yn + k2/2 )
        k4 = F( xn + h   , yn + k3   )
        y(n+1) = 1/6 * h * (k1 + 2*k2 + 2*k3 + k4)
        """
        k1 = np.dot(self.coeff_matrix, self.t_dist)
        k2 = np.dot(self.coeff_matrix, (self.t_dist + 0.5 * self.dt * k1))
        k3 = np.dot(self.coeff_matrix, (self.t_dist + 0.5 * self.dt * k2))
        k4 = np.dot(self.coeff_matrix, (self.t_dist + self.dt * k3))
        self.t_dist = self.t_dist + (k1 + 2 * k2 + 2 * k3 + k4) * self.dt / 6

    def update_pool(self):
        """Update physical properties of pool."""
        if self.mass_water > 0:
            self.vol = (
                self.mass_water * (1 + self.salinity / 1000)
            ) / self.density_water

            if self.angle_wall < (pi / 2):
                # Calculate current depth and surface radius from base radius,
                # wall angle, and volume.
                self.depth = (
                    1
                    / (tan((pi / 2) - self.angle_wall))
                    * (
                        (
                            (
                                3 * self.vol * tan((pi / 2) - self.angle_wall) / pi
                                + self.r_bottom ** 3
                            )
                        )
                        ** (1 / 3)
                        - self.r_bottom
                    )
                )
                self.r_surf = self.r_bottom + self.depth / tan(self.angle_wall)

            elif self.angle_wall == pi / 2:
                self.depth = self.vol / (pi * self.r_bottom ** 2)
                self.r_surf = self.r_bottom
        else:
            # If the pool is dry, set assoc. parameters to empty.
            self.vol = 0
            self.depth = 0
            self.r_surf = self.r_bottom
        self.area_surf = pi * self.r_surf ** 2
