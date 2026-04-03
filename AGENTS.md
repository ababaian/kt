# AGENTS.md

This file defines agent context for `kt` package development and analysis workflow.
If you're an LLM given a task then use the most appropriate 


## Agent Example Structure

```yaml
agent_name:
  description: Brief description of agent purpose
  model: claude-sonnet-4 | claude-opus-4-6 | claude-haiku-4-5  # Optional, defaults to parent
  tools: [tool1, tool2, tool3]  # Available tools for this agent
  isolation: worktree  # Optional, creates isolated git worktree
  context: |
    Additional context or instructions for the agent
```

## Parameter Descriptions

- **description**: 3-5 word summary of agent purpose (required)
- **model**: Specific Claude model to use. Opus for complex reasoning, Sonnet for balanced tasks, Haiku for speed
- **tools**: Array of tools available to agent. Use `[*]` for all tools, or specify: `[Read, Write, Edit, Bash, Grep, Glob]`
- **isolation**: Set to `worktree` to give agent an isolated copy of repo (auto-cleaned if no changes)
- **context**: Additional instructions, domain knowledge, or constraints for the agent

## Example Agents

### Development Agents

```yaml
r-dev:
  description: Develop and validate R package
  tools: [Bash, Read, Grep, Glob]
  context: |
    Focus on R package development best practices:
    - Use devtools/roxygen2 for building/checking
    - Ensure DESCRIPTION dependencies are correct
    - Run R CMD check before suggesting changes
    - Follow Bioconductor guidelines for bioinformatics packages

    R function organization style:
    - Separate each function into its own file named after the main function
    - Group related functions working on the same object in one file
    - File names should be descriptive: load_kt_fasta.R, load_kt_uclust.R
    - Exception: helper functions (%||%, internal utilities) can stay with their main function
    
    KT Network Analysis Design Decisions:
    - create_kt_graph() creates UNDIRECTED graphs (protein similarity)
    - Automatically detects and deduplicates reciprocal alignments (A->B, B->A)
    - Selects best alignment per protein pair (highest bitscore, then identity)
    - Node attributes: family, organism (biological properties)
    - Edge attributes: identity, length, evalue, bitscore (alignment quality)
    
    KT Clustering Analysis Design:
    - Cluster size = number of "H" (hit) records per cluster_id
    - Singletons = clusters with 0 members (only centroid "C" record)
    - Provides log-scale visualization for wide size distributions
    - Family-wise analysis shows clustering patterns across KT families
    - Identity vs size plots reveal clustering quality relationships

r-tester:
  description: Execute R package tests
  tools: [Bash, Read, Write]
  context: |
    Run testthat tests and validate package functionality.
    Focus on edge cases and data integrity for bioinformatics workflows.

r-writer:
  description: Generate R documentation and vignettes
  model: claude-opus-4-6
  tools: [Read, Write, Edit, Grep]
  context: |
    Create roxygen2 documentation and R Markdown vignettes.
    Use bioinformatics terminology appropriately.
    Include working code examples with sample data.
```

### Analysis Agents

```yaml
protein-analyzer:
  description: Analyze KT protein sequences
  model: claude-opus-4-6
  tools: [Read, Write, Bash, Grep]
  context: |
    Expert in protein sequence analysis and bioinformatics.
    Focus on KT protein biodiversity patterns, evolutionary analysis.
    Suggest appropriate R/Bioconductor packages for analysis.

data-validator:
  description: Validate bioinformatics data integrity
  tools: [Read, Bash, Grep]
  context: |
    Validate protein sequence files (FASTA, TSV formats).
    Check for data consistency, missing values, format compliance.
    Ensure sequences are properly formatted and annotated.

visualization-agent:
  description: Create scientific plots and figures
  model: claude-sonnet-4
  tools: [Read, Write, Edit]
  context: |
    Generate publication-ready visualizations using ggplot2.
    Follow scientific visualization best practices.
```

## Usage

Agents are invoked using the Agent tool:
```
Agent(subagent_type="r-package-builder", description="Build package", prompt="Run R CMD check")
```

## Notes

- Agents inherit context and memory from parent conversation
- Use `isolation: worktree` for experimental changes
- Combine multiple agents for complex workflows
- Agents automatically handle R package development conventions