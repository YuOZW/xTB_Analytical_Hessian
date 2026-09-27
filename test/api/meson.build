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

test(
  'C-API',
  executable(
    '@0@-capi-tester'.format(meson.project_name()),
    sources: files('main.c'),
    dependencies: tblite_dep,
  ),
  suite: 'api',
)
# Standalone smoke test for the public GFN1-xTB analytic-Hessian C API.
test(
  'C-API-Hessian',
  executable(
    '@0@-hessian-capi-tester'.format(meson.project_name()),
    sources: files('hessian.c'),
    dependencies: tblite_dep,
  ),
  suite: 'api',
)
