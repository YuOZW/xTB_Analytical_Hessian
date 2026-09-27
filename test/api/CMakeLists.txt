# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with tblite.  If not, see <https://www.gnu.org/licenses/>.

add_executable(
  "${PROJECT_NAME}-api-tester"
  "main.c"
)
target_link_libraries(
  "${PROJECT_NAME}-api-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
)
add_test("C-API" "${PROJECT_NAME}-api-tester")
# Standalone smoke test for the public GFN1-xTB analytic-Hessian C API.
# Kept separate from the broad C-API test so Hessian validation does not
# depend on optional ddX support.
add_executable(
  "${PROJECT_NAME}-hessian-api-tester"
  "hessian.c"
)
target_link_libraries(
  "${PROJECT_NAME}-hessian-api-tester"
  PRIVATE
  "${PROJECT_NAME}-lib"
)
add_test("C-API-Hessian" "${PROJECT_NAME}-hessian-api-tester")
