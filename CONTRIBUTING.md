# Contributing to GenoDot

Thank you for your interest in contributing to GenoDot! This document provides guidelines and information for contributors.

## 🤝 How to Contribute

### Reporting Bugs

- Use the [GitHub Issues](https://github.com/genomics/genodot/issues) page
- Provide a clear, descriptive title
- Include:
  - R version and operating system
  - GenoDot version
  - Sample PAF file (if possible)
  - Command used
  - Error message
  - Expected vs actual behavior

### Suggesting Features

- Open an issue with "Feature Request" label
- Describe the use case and benefits
- Provide examples if possible
- Consider implementation complexity

### Code Contributions

#### 1. Fork the Repository

```bash
# Fork on GitHub, then clone your fork
git clone https://github.com/YOUR_USERNAME/genodot.git
cd genodot
git remote add upstream https://github.com/genomics/genodot.git
```

#### 2. Create a Branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/your-bug-fix
```

#### 3. Make Changes

- Follow the existing code style
- Add comments for complex logic
- Update documentation as needed
- Add tests for new features

#### 4. Test Your Changes

```bash
# Test with sample data
./genodot.R examples/sample.paf

# Run R syntax check
R -e "source('genodot.R')"
```

#### 5. Submit Pull Request

- Push to your fork
- Create pull request on GitHub
- Describe changes clearly
- Link relevant issues
- Wait for review

## 📝 Code Style Guidelines

### R Code Style

- Use 2-space indentation
- Use snake_case for variable names
- Use descriptive function names
- Add Roxygen2 comments for functions
- Keep lines under 80 characters when possible

### Example Function

```r
#' Calculate alignment statistics
#'
#' @param alignments Data frame of alignments
#' @return List with statistics
calc_alignment_stats <- function(alignments) {
  total_alignments <- nrow(alignments)
  total_bp_aligned <- sum(alignments$lenAln)
  avg_identity <- mean(alignments$percentID) * 100
  
  list(
    total_alignments = total_alignments,
    total_bp_aligned = total_bp_aligned,
    avg_identity = avg_identity
  )
}
```

## 🧪 Testing

### Running Tests

```bash
# Basic functionality test
./genodot.R -v

# Test with different options
./genodot.R -p 10 -F png examples/sample.paf
./genodot.R -C Viridis -S examples/sample.paf
```

### Adding Tests

- Create test PAF files in `tests/` directory
- Test edge cases (empty files, large genomes)
- Verify output formats
- Test new features thoroughly

## 📚 Documentation

### Updating README

- Keep installation instructions current
- Update feature list
- Add new examples
- Update badges and links

### Code Documentation

- Use Roxygen2 style comments
- Document all parameters
- Provide return value descriptions
- Include usage examples

## 🔄 Development Workflow

### Before Starting

1. Check existing issues and pull requests
2. Discuss major changes in an issue first
3. Ensure your idea aligns with project goals

### During Development

1. Work in small, focused commits
2. Test frequently
3. Keep code clean and readable
4. Update documentation as you go

### Before Submitting

1. Test all functionality
2. Check code style
3. Update documentation
4. Ensure all tests pass

## 🏷️ Issue Labels

- `bug`: Bug reports and fixes
- `enhancement`: Feature improvements
- `feature`: New features
- `documentation`: Documentation updates
- `good first issue`: Good for new contributors
- `help wanted`: Community help needed
- `priority: high`: High priority issues
- `priority: low`: Low priority issues

## 🎯 Priority Areas

We're currently looking for contributions in:

1. **Performance optimization** for large genomes
2. **Additional color palettes** and themes
3. **Interactive features** and web interface
4. **Integration with other tools** (MUMmer, SyRI, etc.)
5. **Python version** for cross-platform compatibility
6. **Docker containerization**
7. **Galaxy tool wrapper**

## 📋 Release Process

### Version Bumping

- Update version in `genodot.R`
- Update version in README.md
- Create release notes
- Tag release in Git

### Release Checklist

- [ ] All tests pass
- [ ] Documentation updated
- [ ] Version numbers updated
- [ ] Changelog updated
- [ ] GitHub release created
- [ ] Docker image updated (if applicable)

## 🏆 Recognition

Contributors will be:
- Listed in README.md
- Added to AUTHORS file
- Mentioned in release notes
- Invited to join the core team (for significant contributions)

## 📞 Getting Help

- 🐛 [Issues](https://github.com/genomics/genodot/issues)
- 💬 [Discussions](https://github.com/genomics/genodot/discussions)
- 📧 maintainers@genodot.org

## 📄 License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to GenoDot! Your help makes genome visualization better for everyone. 🧬✨
