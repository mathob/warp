#!/usr/bin/env bash

#
# Test script for WholeGenomeGermlineSingleSampleFastq Nextflow pipeline
# This script performs basic validation of the pipeline structure and syntax
#

set -euo pipefail

# Add /Users/mathob/bin to PATH if it exists
if [ -d "/Users/mathob/bin" ]; then
    export PATH="/Users/mathob/bin:$PATH"
fi

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo -e "${GREEN}=== WholeGenomeGermlineSingleSampleFastq Nextflow Pipeline Test ===${NC}"
echo ""

# Function to print status
print_status() {
    local status=$1
    local message=$2
    
    if [ "$status" == "PASS" ]; then
        echo -e "${GREEN}[PASS]${NC} $message"
    elif [ "$status" == "FAIL" ]; then
        echo -e "${RED}[FAIL]${NC} $message"
        exit 1
    elif [ "$status" == "WARN" ]; then
        echo -e "${YELLOW}[WARN]${NC} $message"
    else
        echo -e "[INFO] $message"
    fi
}

# Check if Nextflow is installed
echo "Checking prerequisites..."
NEXTFLOW_FOUND=false
if command -v nextflow >/dev/null 2>&1; then
    NEXTFLOW_VERSION=$(nextflow -version 2>/dev/null | grep version | awk '{print $2}' || echo "unknown")
    print_status "PASS" "Nextflow found (version: $NEXTFLOW_VERSION)"
    NEXTFLOW_FOUND=true
else
    print_status "WARN" "Nextflow not found. Skipping Nextflow-specific tests."
    print_status "INFO" "To install Nextflow: curl -s https://get.nextflow.io | bash"
    NEXTFLOW_VERSION="not_installed"
fi

# Check Nextflow version only if Nextflow is found
if [ "$NEXTFLOW_FOUND" = true ]; then
    REQUIRED_VERSION="21.10.3"
    if [ "$NEXTFLOW_VERSION" != "unknown" ] && [ "$(printf '%s\n' "$REQUIRED_VERSION" "$NEXTFLOW_VERSION" | sort -V | head -n1)" = "$REQUIRED_VERSION" ]; then
        print_status "PASS" "Nextflow version meets minimum requirement ($REQUIRED_VERSION)"
    else
        print_status "WARN" "Nextflow version may be too old or unknown. Required: $REQUIRED_VERSION, Found: $NEXTFLOW_VERSION"
    fi
fi

echo ""
echo "Checking pipeline files..."

# Check main pipeline file
MAIN_PIPELINE="$SCRIPT_DIR/WholeGenomeGermlineSingleSampleFastq.nf"
if [ -f "$MAIN_PIPELINE" ]; then
    print_status "PASS" "Main pipeline file found"
else
    print_status "FAIL" "Main pipeline file not found: $MAIN_PIPELINE"
fi

# Check configuration file
CONFIG_FILE="$SCRIPT_DIR/nextflow.config"
if [ -f "$CONFIG_FILE" ]; then
    print_status "PASS" "Configuration file found"
else
    print_status "FAIL" "Configuration file not found: $CONFIG_FILE"
fi

# Check example parameters file
PARAMS_FILE="$SCRIPT_DIR/params.config"
if [ -f "$PARAMS_FILE" ]; then
    print_status "PASS" "Example parameters file found"
else
    print_status "FAIL" "Example parameters file not found: $PARAMS_FILE"
fi

# Check modules directory
MODULES_DIR="$SCRIPT_DIR/modules"
if [ -d "$MODULES_DIR" ]; then
    print_status "PASS" "Modules directory found"
    
    # Check individual module files
    for module in alignment.nf processing.nf gatk.nf qc.nf; do
        if [ -f "$MODULES_DIR/$module" ]; then
            print_status "PASS" "Module found: $module"
        else
            print_status "FAIL" "Module not found: $module"
        fi
    done
else
    print_status "FAIL" "Modules directory not found: $MODULES_DIR"
fi

echo ""
echo "Validating pipeline syntax..."

# Validate main pipeline syntax only if Nextflow is available
if [ "$NEXTFLOW_FOUND" = true ]; then
    if nextflow config "$MAIN_PIPELINE" >/dev/null 2>&1; then
        print_status "PASS" "Main pipeline syntax is valid"
    else
        print_status "FAIL" "Main pipeline syntax error detected"
    fi
else
    print_status "WARN" "Skipping syntax validation (Nextflow not available)"
fi

# Check for required parameters in main pipeline
echo ""
echo "Checking pipeline structure..."

# Check if main pipeline has required elements
if grep -q "nextflow.enable.dsl.*2" "$MAIN_PIPELINE"; then
    print_status "PASS" "DSL2 syntax enabled"
else
    print_status "WARN" "DSL2 syntax not explicitly enabled"
fi

if grep -q "workflow {" "$MAIN_PIPELINE"; then
    print_status "PASS" "Main workflow block found"
else
    print_status "FAIL" "Main workflow block not found"
fi

if grep -q "include.*from.*modules" "$MAIN_PIPELINE"; then
    print_status "PASS" "Module imports found"
else
    print_status "FAIL" "Module imports not found"
fi

# Check configuration file syntax
echo ""
echo "Validating configuration..."

if [ -f "$CONFIG_FILE" ]; then
    # Simple syntax check for Groovy/Nextflow config
    if grep -q "params {" "$CONFIG_FILE" && grep -q "process {" "$CONFIG_FILE"; then
        print_status "PASS" "Configuration file structure looks correct"
    else
        print_status "WARN" "Configuration file may have structural issues"
    fi
    
    # Check for profiles
    if grep -q "profiles {" "$CONFIG_FILE"; then
        print_status "PASS" "Execution profiles defined"
    else
        print_status "WARN" "No execution profiles found in config"
    fi
fi

echo ""
echo "Checking container compatibility..."

# Check Docker availability
if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
    print_status "PASS" "Docker is available"
else
    print_status "WARN" "Docker not available or not running"
fi

# Check Singularity availability  
if command -v singularity >/dev/null 2>&1 || command -v apptainer >/dev/null 2>&1; then
    print_status "PASS" "Singularity/Apptainer is available"
else
    print_status "WARN" "Singularity/Apptainer not available"
fi

# Check Conda availability
if command -v conda >/dev/null 2>&1 || command -v mamba >/dev/null 2>&1; then
    print_status "PASS" "Conda/Mamba is available"
else
    print_status "WARN" "Conda/Mamba not available"
fi

echo ""
echo "Performing dry run test..."

if [ "$NEXTFLOW_FOUND" = true ]; then
    # Test syntax validation using nextflow config command
    if nextflow config "$MAIN_PIPELINE" >/dev/null 2>&1; then
        print_status "PASS" "Pipeline configuration validation successful"
    else
        print_status "WARN" "Pipeline configuration validation encountered issues"
    fi
    
    # Test help option to validate parameter parsing
    if nextflow run "$MAIN_PIPELINE" --help >/dev/null 2>&1; then
        print_status "PASS" "Pipeline help command works correctly"
    else
        print_status "WARN" "Pipeline help command encountered issues"
    fi
else
    print_status "WARN" "Skipping dry run test (Nextflow not available)"
fi

echo ""
echo -e "${GREEN}=== Test Summary ===${NC}"
print_status "INFO" "Pipeline structure validation completed"
print_status "INFO" "Review any warnings above before running with real data"

echo ""
echo "Next steps:"
echo "1. Edit params.config with your actual file paths"
echo "2. Choose an appropriate execution profile (docker, singularity, conda)"
echo "3. Run the pipeline: nextflow run WholeGenomeGermlineSingleSampleFastq.nf -c params.config -profile <profile>"

echo ""
print_status "INFO" "Test completed successfully!"
