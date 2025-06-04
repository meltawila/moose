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

#include "ThermalTungstenProperties.h"
#include "libmesh/utility.h"

registerMooseObject("SolidPropertiesApp", ThermalTungstenProperties);

InputParameters
ThermalTungstenProperties::validParams()
{
  InputParameters params = ThermalSolidProperties::validParams();
  params.addClassDescription("Tungsten thermal properties.");
  return params;
}

ThermalTungstenProperties::ThermalTungstenProperties(const InputParameters & parameters)
  : ThermalSolidProperties(parameters),
   _kN(8.7401e-01),
   _kA0(7.348e+05),
   _kA1(2.544e+01),
   _kA2(-8.304e+03),
   _kA3(1.180e+06),
   _kB0(-3.679),
   _kB1(1.181e+02),
   _kB2(5.879e+01),
   _kB3(2.867),
   _kC0(-2.052e-02),
   _kC1(4.741e-01),
   _cN(3.030),
   _cA0(3.103e+02),
   _cA1(-8.815),
   _cA2(1.295e+02),
   _cA3(1.874e+03),
   _cB0(1.301e-01),
   _cB1(2.225e-02),
   _cB2(-7.224e-03),
   _cB3(3.539e-03),
   _cB_2(-3.061e-04)
{
}

Real
ThermalTungstenProperties::cp_from_T(const Real & T) const
{
  if (T < 293)
  return _cA0 * Utility::pow<_cN>(T / 1000) / (1 + _cA1 * T / 1000 + _cA2 * Utility::pow<2>(T / 1000) + _cA3 * Utility::pow<3>(T / 1000));
  else
  return _cB0 + _cB1 * T / 1000 + _cB2 * Utility::pow<2>(T / 1000) + _cB3 * Utility::pow<3>(T / 1000) + _cB_2 / Utility::pow<2>(T / 1000);
}

void
ThermalTungstenProperties::cp_from_T(const Real & T, Real & cp, Real & dcp_dT) const
{
  cp = cp_from_T(T);
  if (T < 293) {
    dcp_dT = 1 / 1000 * (_cA0 * _cN * Utility::pow<_cN - 1>(T / 1000) * (1 + _cA1 * T / 1000 + _cA2 * Utility::pow<2>(T / 1000) + _cA3 * Utility::pow<3>(T / 1000)) 
            - _cA0 * Utility::pow<_cN>(T / 1000) * (_cA1 + 2 * _cA2 * T / 1000 + 3 * _cA3 * Utility::pow<2>(T / 1000)))
            / Utility::pow<2>(1+ _cA1 * T / 1000 + _cA2 * Utility::pow<2>(T / 1000) + _cA3 Utility::pow<3>(T / 1000));}
  else {
    dcp_dT = 1 / 1000 * ( _cB1 + 2 * _cB2 * T / 1000 + 3 * _cB3 * Utility::pow<2>(T / 1000) - 2 * _cB_2 / Utility::pow<3>(T / 1000));}
}

Real
ThermalTungstenProperties::k_from_T(const Real & T) const
{
  if (T < 55)
  return _kA0 * Utility::pow<_kN>(T / 1000) / (1 + _kA1 * T / 1000 + _kA2 * Utility::pow<2>(T / 1000) + _kA3 * Utility::pow<3>(T / 1000));
  else
  return (_kB0 + _kB1 * T / 1000 + _kB2 * Utility::pow<2>(T / 1000) + _kB3 * Utility::pow<3>(T / 1000))/(_kC0 + _kC1 * T / 1000+ Utility::pow<2>(T / 1000));
}

void
ThermalTungstenProperties::k_from_T(const Real & T, Real & k, Real & dk_dT) const
{
  k = k_from_T(T);
  if (T < 55) {
    dk_dT = 1/1000 * (_kA0 * _kN * Utility::pow<_kN-1>(T / 1000) * (1 + _kA1 * T / 1000 + _kA2 * Utility::pow<2>(T / 1000) + _kA3 * Utility::pow<3>(T / 1000)) 
            - _kA0 * Utility::pow<_kN>(T / 1000)*(_kA1 + 2 * _kA2 * T / 1000 + 3 * _kA3 * Utility::pow<2>(T / 1000)))
            / Utility::pow<2>(1+ _kA1 * T / 1000 + _kA2 * Utility::pow<2>(T / 1000) + _kA3 Utility::pow<3>(T / 1000))}
  else {
    dk_dT = 1/1000 * ((_kB1 + 2 * _kB2 * T / 1000 + 3 * _kB3 * Utility::pow<2>(T / 1000)) * (_kC0 + _kC1 * T / 1000 + Utility::pow<2>(T / 1000))
            - (_kB0 + _kB1 * T / 1000 + _kB2 * Utility::pow<2>(T / 1000)+ _kB3 * Utility::pow<3>(T / 1000)) * (_kC1 + 2 * T / 1000))
            / Utility::pow<2>(_kC0+ _kC1 * T / 1000 + _kC2 * Utility::pow<2>(T / 1000))}
}

Real
ThermalTungstenProperties::rho_from_T(const Real & T) const
{
  if (T <= 294)
    return 19250 / Utility::pow<3>(1 + (-8.529e-02 - 9.915e-02 * T / 1000 + 2.257 * Utility::pow<2>(T / 1000) -3.157 * Utility::pow<3>(T / 1000)) / 100);
  else
    return 19250 / Utility::pow<3>(1 + (-1.400e-01 + 4.869e-01 * T / 1000 - 3.056e-02 * Utility::pow<2>(T / 1000) + 2.234e-02 * Utility::pow<3>(T / 1000)) / 100);
}

void
ThermalTungstenProperties::rho_from_T(const Real & T, Real & rho, Real & drho_dT) const
{
  rho = rho_from_T(T);
  if (T <= 294)
    drho_dT = -3 * 19250e3 * (-9.915e-02 + 2 * 2.257 * T / 1000 + 3 * -3.157 * Utility::pow<2>(T / 1000)) 
              / Utility::pow<4>(100 - 8.529e-02 - 9.915e-02 * T / 1000 + 2.257 * Utility::pow<2>(T / 1000) - 3.157 * Utility::pow<3>(T / 1000));
  else
    drho_dT = -3 * 19250e3 * (4.869e-01 + 2 * -3.056e-02 * T / 1000 + 3 * 2.234e-02 * Utility::pow<2>(T / 1000)) 
              / Utility::pow<4>(100 - 1.400e-01 + 4.869e-01 * T / 1000 - 3.056e-02 * Utility::pow<2>(T / 1000) + 2.234e-02 * Utility::pow<3>(T / 1000));
}
