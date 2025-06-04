//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

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
  const Real _kN;
  const Real _kA0;
  const Real _kA1;
  const Real _kA2;
  const Real _kA3;
  const Real _kB0;
  const Real _kB1;
  const Real _kB2;
  const Real _kB3;
  const Real _kC0;
  const Real _kC1;
  const Real _cN;
  const Real _cA0;
  const Real _cA1;
  const Real _cA2;
  const Real _cA3;
  const Real _cB0;
  const Real _cB1;
  const Real _cB2;
  const Real _cB3;
  const Real _cB_2;
};

#pragma GCC diagnostic pop
