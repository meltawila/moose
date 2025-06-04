//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Data from: Milner, Justin L., Peter Karkos, and Jessica Jane Bowers. 
//* Space Nuclear Propulsion (SNP) Material Property Handbook. 
//* No. SNP-HDBK-0008. National Aeronautics and Space Administration, 2024.

#pragma once

#include "ThermalSolidProperties.h"

/**
 * Tungsten thermal solid properties as a function of temperature.
 */
class ThermalTungstenProperties : public ThermalSolidProperties
{
public:
  static InputParameters validParams();

  ThermalTungstenProperties(const InputParameters & parameters);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

  virtual Real k_from_T(const Real & T) const override;

  virtual void k_from_T(const Real & T, Real & k, Real & dk_dT) const override;

  virtual Real cp_from_T(const Real & T) const override;

  virtual void cp_from_T(const Real & T, Real & cp, Real & dcp_dT) const override;

  virtual Real rho_from_T(const Real & T) const override;

  virtual void rho_from_T(const Real & T, Real & rho, Real & drho_dT) const override;

protected:
  // Constants used in specific heat and thermal conductivity relations
  static constexpr Real _kN  = 8.7401e-01;
  static constexpr Real _kA0  = 7.348e+05;
  static constexpr Real _kA1  = 2.544e+01;
  static constexpr Real _kA2  = -8.304e+03;
  static constexpr Real _kA3  = 1.180e+06;
  static constexpr Real _kB0  = -3.679;
  static constexpr Real _kB1  = 1.181e+02;
  static constexpr Real _kB2  = 5.879e+01;
  static constexpr Real _kB3  = 2.867;
  static constexpr Real _kC0  = -2.052e-02;
  static constexpr Real _kC1  = 4.741e-01;
  static constexpr Real _cN   = 3.030;
  static constexpr Real _cA0  = 3.103e+02;
  static constexpr Real _cA1  = -8.815;
  static constexpr Real _cA2  = 1.295e+02;
  static constexpr Real _cA3  = 1.874e+03;
  static constexpr Real _cB0  = 1.301e-01;
  static constexpr Real _cB1  = 2.225e-02;
  static constexpr Real _cB2  = -7.224e-03;
  static constexpr Real _cB3  = 3.539e-03;
  static constexpr Real _cB_2  = -3.061e-04;
};

#pragma GCC diagnostic pop
