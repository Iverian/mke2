#ifndef MKE2_INCLUDE_EXPORT_HPP_
#define MKE2_INCLUDE_EXPORT_HPP_

#include "triangulation.hpp"
#include "vec.hpp"

#include <ostream>

void mv2_export(std::ostream& os, const Triangulation& t, const Vec& values);
void csv_export(std::ostream& os, const Triangulation& t, const Vec& values);

#endif // MKE2_INCLUDE_MV2_EXPORT_HPP_