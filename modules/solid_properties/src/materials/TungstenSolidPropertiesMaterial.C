//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TungstenSolidPropertiesMaterial.h"
#include "ThermalSolidProperties.h"
#include "SolidPropertiesNames.h"

registerMooseObject("SolidPropertiesApp", TungstenSolidPropertiesMaterial);
registerMooseObject("SolidPropertiesApp", ADTungstenSolidPropertiesMaterial);

template <bool is_ad>
InputParameters
TungstenSolidPropertiesMaterialTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();
  params.addRequiredCoupledVar("temperature", "Temperature");
  params.addRequiredParam<UserObjectName>("sp", "The name of the user object for solid properties");
  params.addParam<std::string>(SolidPropertiesNames::specific_heat,
                               SolidPropertiesNames::specific_heat,
                               "Name to be used for the isobaric specific heat");
  params.addParam<std::string>(SolidPropertiesNames::thermal_conductivity,
                               SolidPropertiesNames::thermal_conductivity,
                               "Name to be used for the thermal conductivity");
  params.addParam<std::string>(SolidPropertiesNames::density,
                               SolidPropertiesNames::density,
                               "Name to be used for the density");
  params.addClassDescription("Computes solid thermal properties as a function of temperature");
  return params;
}

template <bool is_ad>
TungstenSolidPropertiesMaterialTempl<is_ad>::TungstenSolidPropertiesMaterialTempl(
    const InputParameters & parameters)
  : Material(parameters),
    _temperature(coupledGenericValue<is_ad>("temperature")),

    _cp(declareGenericProperty<Real, is_ad>(
        getParam<std::string>(SolidPropertiesNames::specific_heat))),
    _k(declareGenericProperty<Real, is_ad>(
        getParam<std::string>(SolidPropertiesNames::thermal_conductivity))),
    _rho(declareGenericProperty<Real, is_ad>(getParam<std::string>(SolidPropertiesNames::density))),

    _Tungsten(getUserObject<TungstenProperties>("sp"))
{
}

template <bool is_ad>
void
TungstenSolidPropertiesMaterialTempl<is_ad>::computeQpProperties()
{
  _cp[_qp] = _sp.cp_from_T(_temperature[_qp]);
  _k[_qp] = _sp.k_from_T(_temperature[_qp]);
  _rho[_qp] = _sp.rho_from_T(_temperature[_qp]);
}

template class TungstenSolidPropertiesMaterialTempl<false>;
template class TungstenSolidPropertiesMaterialTempl<true>;
