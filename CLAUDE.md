# IntegrateALL Pipeline - Claude Code Configuration

## Project Overview
IntegrateALL is an advanced, high-performance B-cell acute lymphoblastic leukemia (B-ALL) RNA-seq analysis pipeline. This is an improved version of the original "Blast-o-Matic-Fusioninator_cluster" pipeline with unified Python codebase, modular architecture, and enhanced performance.

## Key Technologies
- **Snakemake**: Workflow management system
- **Python 3.9+**: Primary programming language
- **R 4.0+**: For specific bioinformatics tools (ALLCatchR, RNASeqCNV)
- **Conda/Mamba**: Environment management
- **STAR**: RNA-seq alignment
- **Arriba/FusionCatcher**: Fusion detection
- **GATK**: Variant calling
- **ALLCatchR**: B-ALL subtype classification

## Project Structure

```
IntegrateALL_pipeline/
├── workflow/
│   ├── Snakefile           # Main workflow entry point
│   ├── rules/              # Modular Snakemake rules
│   │   ├── common.smk      # Common functions and utilities
│   │   ├── qc.smk          # Quality control (FastQC, MultiQC)
│   │   ├── alignment.smk   # STAR alignment
│   │   ├── fusion_detection.smk  # Arriba & FusionCatcher
│   │   ├── variant_calling.smk   # GATK pipeline
│   │   ├── cnv_analysis.smk      # RNASeqCNV analysis
│   │   ├── classification.smk    # ALLCatchR & integration
│   │   └── reporting.smk         # Report generation
│   ├── scripts/            # Python analysis scripts
│   │   ├── integrated_classification.py  # Main classification logic
│   │   ├── integrate_fusion_results.py   # Fusion integration
│   │   ├── generate_final_report.py      # HTML report generation
│   │   └── ...
│   └── envs/              # Conda environment files
├── config/
│   ├── config.yaml        # Main pipeline configuration
│   ├── samples.tsv        # Sample metadata
│   └── resources.yaml     # Resource specifications
├── resources/             # Reference data and databases
├── results/               # Pipeline outputs
└── docs/                  # Documentation
```

## Development Environment

### Python Environment
- **Python Version**: 3.9+
- **Main Dependencies**:
  - pandas >= 1.5.0
  - numpy >= 1.21.0
  - scikit-learn >= 1.1.0
  - matplotlib >= 3.5.0
  - seaborn >= 0.11.0
  - pyyaml >= 6.0
  - snakemake >= 7.0

### IDE Configuration
- **Recommended IDE**: VS Code with Python extension
- **Linting**: flake8, black (code formatting)
- **Type Checking**: mypy
- **Testing**: pytest

## Development Commands

### Environment Setup
```bash
# Create and activate main environment
conda env create -f workflow/envs/python.yaml
conda activate integrateall_python

# Install additional development tools
pip install pytest black flake8 mypy
```

### Running the Pipeline
```bash
# Local execution
snakemake --cores 16 --use-conda

# Dry run (check workflow)
snakemake --dry-run --cores 1

# Run specific target
snakemake results/reports/sample_001/sample_001_final_report.html --cores 8

# Cluster execution (SLURM)
snakemake --profile config/slurm --jobs 50
```

### Development and Testing
```bash
# Run Python tests
pytest .tests/

# Format Python code
black workflow/scripts/

# Lint Python code
flake8 workflow/scripts/

# Type checking
mypy workflow/scripts/

# Profile performance
snakemake --benchmark-extended --cores 8
```

## Key Python Modules

### `integrated_classification.py`
**Purpose**: Combines results from multiple analysis tools (ALLCatchR, karyotype prediction, fusion detection, CNV analysis) for comprehensive B-ALL subtype classification.

**Key Classes**:
- `IntegratedClassifier`: Main classification logic
- Methods for loading and integrating multi-modal data

**Usage Example**:
```python
classifier = IntegratedClassifier(rules_file, confidence_threshold=0.7)
classification, confidence, evidence = classifier.classify_sample(
    allcatchr_results, karyotype_results, fusion_results, 
    cnv_results, hotspot_results
)
```

### `integrate_fusion_results.py`
**Purpose**: Integrates and filters fusion calls from Arriba and FusionCatcher.

**Key Classes**:
- `FusionIntegrator`: Fusion integration and filtering logic
- Methods for quality filtering and driver gene identification

### `generate_final_report.py`  
**Purpose**: Creates comprehensive HTML reports with interactive visualizations.

**Key Classes**:
- `ReportGenerator`: HTML report generation with matplotlib plots
- Methods for data visualization and report formatting

## Configuration Files

### `config/config.yaml`
Main pipeline configuration including:
- Reference genome paths
- Tool parameters
- Resource allocations
- Database locations

### `config/samples.tsv`
Sample metadata in TSV format:
```tsv
sample_id	fastq1	fastq2	condition	batch
sample_001	/path/to/R1.fastq.gz	/path/to/R2.fastq.gz	B-ALL	batch1
```

## Workflow Rules

### Core Analysis Steps
1. **Quality Control** (`qc.smk`): FastQC and MultiQC analysis
2. **Alignment** (`alignment.smk`): STAR RNA-seq alignment
3. **Fusion Detection** (`fusion_detection.smk`): Arriba and FusionCatcher
4. **Variant Calling** (`variant_calling.smk`): GATK best practices
5. **CNV Analysis** (`cnv_analysis.smk`): RNASeqCNV chromosomal analysis
6. **Classification** (`classification.smk`): ALLCatchR and integration
7. **Reporting** (`reporting.smk`): Final HTML report generation

### Resource Management
- Dynamic resource allocation based on rule requirements
- Cluster-specific configurations for SLURM/PBS
- Automatic temporary file cleanup
- Checkpoint-based resumable execution

## Testing

### Unit Tests
Located in `.tests/` directory:
- Test individual Python functions
- Mock data for reproducible testing
- Integration tests for complete workflows

### Validation
- Test datasets with known classifications
- Benchmark against original pipeline results
- Performance regression testing

## Debugging and Troubleshooting

### Common Issues
1. **Memory Issues**: Adjust memory settings in `config.yaml`
2. **Conda Environment**: Clear cache and recreate environments
3. **Reference Files**: Verify paths and permissions
4. **Cluster Configuration**: Check SLURM settings and quotas

### Logging
- Pipeline logs in `logs/` directory
- Rule-specific log files
- Snakemake execution logs
- Benchmark performance data

### Development Debugging
```bash
# Verbose Snakemake output
snakemake --debug-dag --verbose

# Python script debugging
python -m pdb workflow/scripts/script_name.py

# Interactive Python debugging
import pdb; pdb.set_trace()
```

## Performance Considerations

### Optimization Tips
- Use SSD storage for temporary files
- Optimize thread allocation per rule
- Monitor memory usage with system tools
- Use cluster scheduling for large datasets

### Benchmarking
- Built-in Snakemake benchmarking
- Resource usage tracking
- Performance comparison with original pipeline
- Scalability testing with multiple samples

## Contributing

### Code Style
- Follow PEP 8 for Python code
- Use type hints for function signatures
- Comprehensive docstrings for all functions
- Consistent error handling patterns

### Pull Request Process
1. Create feature branch
2. Add tests for new functionality
3. Update documentation
4. Run full test suite
5. Submit pull request with detailed description

## Deployment

### Production Setup
- Use dedicated compute resources
- Configure cluster-specific settings
- Set up monitoring and alerting
- Regular reference database updates

### Container Deployment
- Singularity containers supported
- Docker images available
- Reproducible execution environments

This pipeline represents a significant advancement in B-ALL analysis workflows, providing robust, scalable, and maintainable bioinformatics analysis capabilities.