"""Helper functions for running Snakemake."""

__author__ = "Kat Steinke"

#  Copyright (c) 2026 Kat Steinke
#     This program is distributed under version 3 of the GNU General Public License.
#      You should have received a copy of the GNU General Public License
#        along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

import pathlib


def is_gzipped(fastq_dir: pathlib.Path, barcode_number: str) -> bool:
    """Determine whether all files in the barcode directory are gzipped.

    Arguments:
        fastq_dir:      the directory containing barcodeXX subdirectories for the run
        barcode_number: the barcode number for which fastq file format should be identified
                        (as one of either .fastq or .fastq.gz)

    Returns:
        True if fastq files are gzipped, False otherwise
    Raises:
        FileNotFoundError:  if the supplied fastq directory does not contain the required barcode
                            directory
        ValueError:         if the extension is unsupported or there are multiple formats in the
                            directory
    """
    barcode_dir = fastq_dir / f"barcode{barcode_number}"
    if not barcode_dir.exists():
        raise FileNotFoundError(f"Barcode directory barcode{barcode_number} "
                                f"not found in {fastq_dir}.")
    extensions = [read_file.suffix for read_file in barcode_dir.iterdir() if read_file.is_file()]
    if all(extension == ".gz" for extension in extensions):
        return True
    if all(extension == ".fastq" for extension in extensions):
        return False
    # if we haven't returned by now something is wrong
    raise ValueError(f"Extensions {sorted(extensions)} not supported.")

