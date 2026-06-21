from file_parser import Parse_section_tsv

#metrics_file = "in/LocalApp2_MetricsOutput.tsv"
metrics_file = "in/Dragen_V2_MetricsOutput.tsv"

headers, sections = Parse_section_tsv(
    metrics_file,
    key_value_sections=[]
)

print("HEADERS")
for header in headers:
    print(header)

print("\nSECTIONS")
for section_name, df in sections.items():
    print(f"\n[{section_name}]")
    print(f"Shape: {df.shape}")
    print("Columns:", df.columns)
    print(df)