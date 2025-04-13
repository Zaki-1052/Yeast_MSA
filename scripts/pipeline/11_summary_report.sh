#!/bin/bash

echo "=== Generating final summary report ==="

# Generate method comparison summary
echo "Method Comparison:" > "results/merged/summary/method_comparison.txt"
for GROUP in WT-37 WTA STC CAS STC-vs-STCCTRL CAS-vs-CASCTRL; do
    if [ -f "results/merged/summary/${GROUP}/highconf_count.txt" ] && [ -f "results/merged/summary/${GROUP}/direct_highconf_count.txt" ]; then
        MERGED_COUNT=$(cat "results/merged/summary/${GROUP}/highconf_count.txt")
        DIRECT_COUNT=$(cat "results/merged/summary/${GROUP}/direct_highconf_count.txt")
        echo "  $GROUP: $MERGED_COUNT high-conf variants (merged), $DIRECT_COUNT high-conf variants (direct)" >> "results/merged/summary/method_comparison.txt"
    else
        echo "  $GROUP: Incomplete data" >> "results/merged/summary/method_comparison.txt"
    fi
done

# Generate the full analysis report
echo "=== Variant Analysis Summary Report ===" > "results/merged/summary/analysis_report.txt"
echo "Date: $(date)" >> "results/merged/summary/analysis_report.txt"
echo "" >> "results/merged/summary/analysis_report.txt"

echo "Treatment-Specific Variants:" >> "results/merged/summary/analysis_report.txt"
for GROUP in WT-37 WTA STC CAS STC-vs-STCCTRL CAS-vs-CASCTRL; do
    if [ -f "results/merged/summary/${GROUP}/specific_count.txt" ] && [ -f "results/merged/summary/${GROUP}/highconf_count.txt" ]; then
        TOTAL=$(cat "results/merged/summary/${GROUP}/specific_count.txt")
        HC=$(cat "results/merged/summary/${GROUP}/highconf_count.txt")
        echo "  $GROUP: $TOTAL total, $HC high-confidence" >> "results/merged/summary/analysis_report.txt"
    else
        echo "  $GROUP: No data available" >> "results/merged/summary/analysis_report.txt"
    fi
done

echo "" >> "results/merged/summary/analysis_report.txt"
echo "Treatment-Unique Variants:" >> "results/merged/summary/analysis_report.txt"
cat "results/merged/summary/unique_variants.txt" | sed 's/^/  /' >> "results/merged/summary/analysis_report.txt"

# Add consistency information to the report
echo "" >> "results/merged/summary/analysis_report.txt"
echo "Variant Consistency Within Replicates:" >> "results/merged/summary/analysis_report.txt"
for GROUP in WT-37 WTA STC CAS; do
    if [ -f "results/merged/summary/${GROUP}/consistency.txt" ]; then
        cat "results/merged/summary/${GROUP}/consistency.txt" | sed 's/^/  /' >> "results/merged/summary/analysis_report.txt"
        echo "" >> "results/merged/summary/analysis_report.txt"
    fi
done

# Add method comparison to the report
echo "Variant Analysis Method Comparison:" >> "results/merged/summary/analysis_report.txt"
cat "results/merged/summary/method_comparison.txt" | sed 's/^/  /' >> "results/merged/summary/analysis_report.txt"

# Add treatment group comparison results
echo "" >> "results/merged/summary/analysis_report.txt"
echo "Treatment Group Comparisons:" >> "results/merged/summary/analysis_report.txt"
for COMPARISON in STC-vs-WTA CAS-vs-WT37 WT37-vs-WTA; do
    if [ -f "results/merged/summary/${COMPARISON}/comparison.txt" ]; then
        echo "  $COMPARISON:" >> "results/merged/summary/analysis_report.txt"
        cat "results/merged/summary/${COMPARISON}/comparison.txt" | sed 's/^/    /' >> "results/merged/summary/analysis_report.txt"
        echo "" >> "results/merged/summary/analysis_report.txt"
    fi
done

echo "Summary report available at: results/merged/summary/analysis_report.txt"