# kt: The KT Protein Universe

`kt` package is an analysis of KT protein biodiversity across the SRA-Logan
Assemblage. This package is a wrapper for protein sequence analysis and 
evolution of KT proteins.

> **⚠️ Development Status**: This package is currently in early development.

## Installation

Development version of `kt`

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install kt from GitHub
devtools::install_github("ababaian/kt")
```

## Quick Start (see: `kt_analysis.Rmd`)

## Data Manifest

- **Source Data**: `kt1.M1.tsv` (28,351 protein sequences)
- **Size**: ~39MB of bioinformatic data
- **Format**: Tab-separated values with protein metadata and sequences
- **Usage**: Data will be processed and made available through package functions

### Data Access

```r
# When implemented:
# Load processed data
# data("kt_data")

# Or load from raw file
# raw_data <- load_kt_data("kt1.M1.tsv")
```

## Documentation

- **Function Reference**: See package documentation (`?kt`)
- **Vignettes**: Comprehensive usage guides (planned)
- **Examples**: All functions include working examples (planned)

## Development

This package is actively under development. To contribute:

1. **Check the TODO**: See [TODO.md](TODO.md) for current development tasks
2. **Report Issues**: Use GitHub issues for bug reports and feature requests
3. **Contribute**: Fork the repository and submit pull requests

### Development Setup

```r
# Clone the repository
git clone https://github.com/username/kt.git

# Install development dependencies
devtools::install_dev_deps()

# Load the package for development
devtools::load_all()

# Run tests
devtools::test()

# Check package
devtools::check()
```

## Requirements

- **R Version**: R >= 4.0.0
- **Dependencies**: Will be listed in DESCRIPTION file as development progresses
- **System**: Cross-platform (Windows, macOS, Linux)

## Author and Maintainer

**Artem Babaian**
Email: a.babaian@utoronto.ca
Affiliation: University of Toronto

## License

This package is licensed under the GNU Affero General Public License v3.0 (AGPL-3).

See [LICENSE](LICENSE) and [LICENSE.md](LICENSE.md) for full license text.

### License Summary

The AGPL-3 license allows you to:
- ✅ Use the software for any purpose
- ✅ Study and modify the source code
- ✅ Share the software and modifications

With these requirements:
- 📋 Include the original license and copyright notice
- 📋 State changes made to the code
- 📋 Make source code available when distributing
- 📋 **Network use triggers copyleft** (AGPL-3 specific)

## Citation

If you use this package in your research, please cite:

```
Babaian, A. (2024). kt: Bioinformatic Analysis of KT Protein Biodiversity.
R package version 0.1.0.
```

(BibTeX and formal citation will be added upon publication)

## Acknowledgments

- Contributors and collaborators (to be added)
- Data sources and collaborating institutions
- Funding sources (if applicable)

## Project Status

- 🔴 **Core Functions**: Not implemented
- 🔴 **Documentation**: Basic structure only
- 🔴 **Tests**: Placeholder only
- 🔴 **Vignettes**: Not created
- 🟡 **Package Structure**: Basic framework complete
- 🟢 **License**: Complete
- 🟢 **Repository**: Initialized

See [TODO.md](TODO.md) for detailed development roadmap.

---

**Last Updated**: March 30, 2026
**Package Version**: 0.1.0
**Development Stage**: Initial Development