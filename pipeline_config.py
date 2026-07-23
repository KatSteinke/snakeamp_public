"""Load default pipeline config."""

#  Copyright (c) 2026 Kat Steinke
#     This program is distributed under version 3 of the GNU General Public License.
#      You should have received a copy of the GNU General Public License
#        along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import pathlib

import yaml

# import check_config - TODO: what config checks do we need in this one?

default_config_file = pathlib.Path(__file__).parent / "config" / "pipeline_routine.yml"

with open(default_config_file, "r", encoding="utf-8") as default_config:
    WORKFLOW_DEFAULT_CONF = yaml.safe_load(default_config)

#check_config.check_config(WORKFLOW_DEFAULT_CONF)
