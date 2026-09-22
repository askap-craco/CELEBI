import os
import sys
import re
import glob

def parse_nf_files(directory):
    """Scans all .nf files to map processes to their labels."""
    process_map = {}
    print(f"Scanning .nf files in: {directory}...")
    
    nf_files = glob.glob(os.path.join(directory, '**/*.nf'), recursive=True)
    if not nf_files:
        print("  -> No .nf files found!")
        return process_map

    for filepath in nf_files:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Split file by the word 'process ' to isolate process blocks
        process_blocks = re.split(r'\bprocess\s+', content)[1:]
        
        for block in process_blocks:
            # Extract the process name
            name_match = re.match(r'([a-zA-Z0-9_]+)', block)
            if not name_match:
                continue
            proc_name = name_match.group(1)
            
            # Extract all labels assigned to this process
            # Matches: label 'foo' OR label "foo"
            labels = re.findall(r'label\s+[\'"]([a-zA-Z0-9_]+)[\'"]', block)
            process_map[proc_name] = labels

    print(f"  -> Found {len(process_map)} processes across {len(nf_files)} files.\n")
    return process_map

def parse_config_file(filepath):
    """Extracts clusterOptions declarations from a config file."""
    config_data = {'withName': {}, 'withLabel': {}}
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"Error: Could not read config file {filepath}")
        return config_data

    current_context = None
    brace_depth = 0

    for line_num, line in enumerate(lines, 1):
        clean_line = line.strip()
        if clean_line.startswith('//'):
            continue

        # Detect entering a withName or withLabel block
        if clean_line.startswith('withName') or clean_line.startswith('withLabel'):
            # Extract everything between the colon and the opening brace
            match = re.search(r'(withName|withLabel)\s*:\s*[\'"]?(.*?)[\'"]?\s*\{', clean_line)
            if match:
                type_ = match.group(1)
                # Handle regex/piped names like "do_ref_correlation|do_correlation"
                targets = [t.strip() for t in match.group(2).split('|')]
                current_context = (type_, targets)

        brace_depth += clean_line.count('{')
        brace_depth -= clean_line.count('}')
        
        # Reset context when we exit the block
        if brace_depth <= 1:
            current_context = None

        # Record clusterOptions if we are inside a context
        if current_context and 'clusterOptions' in clean_line and '=' in clean_line:
            type_, targets = current_context
            for target in targets:
                # Store the line number where this was defined
                # If multiple configs define it, this tracks the last one (which matches Nextflow's behavior)
                config_data[type_][target] = f"{os.path.basename(filepath)}:Line {line_num}"

    return config_data

def lint_pipeline(directory, config_paths):
    # 1. Map all processes and their labels
    process_map = parse_nf_files(directory)
    
    # 2. Parse all provided config files
    combined_config = {'withName': {}, 'withLabel': {}}
    for conf in config_paths:
        data = parse_config_file(conf)
        combined_config['withName'].update(data['withName'])
        combined_config['withLabel'].update(data['withLabel'])

    # 3. Detect Collisions
    print("--- Checking for Configuration Collisions ---")
    collisions_found = False

    for proc_name, labels in process_map.items():
        # Check if the process has a specific withName clusterOptions
        has_name_opt = combined_config['withName'].get(proc_name)
        
        # Check which of its labels define clusterOptions
        label_opts = {}
        for label in labels:
            if label in combined_config['withLabel']:
                label_opts[label] = combined_config['withLabel'][label]

        # Case 1: withName vs. withLabel(s)
        if has_name_opt and label_opts:
            collisions_found = True
            print(f"[COLLISION: Name vs Label] Process '{proc_name}'")
            print(f"  -> Gets clusterOptions from withName ({has_name_opt})")
            for lbl, loc in label_opts.items():
                print(f"  -> AND from withLabel: {lbl} ({loc})")
            print("  ! RESULT: Nextflow uses withName and SILENTLY IGNORES the label(s).\n")

        # Case 2: Multiple Labels colliding
        if len(label_opts) > 1:
            collisions_found = True
            print(f"[COLLISION: Label vs Label] Process '{proc_name}'")
            print(f"  -> Has multiple labels defining clusterOptions:")
            for lbl, loc in label_opts.items():
                print(f"     - '{lbl}' ({loc})")
            print("  ! RESULT: Unpredictable. Nextflow will silently drop all but one of these.\n")

    if not collisions_found:
        print("Success! No withName vs. withLabel or Label vs. Label collisions detected.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python nf_pipeline_linter.py <directory_with_nf_files> <config_file_1> [config_file_2 ...]")
        sys.exit(1)
    
    pipeline_dir = sys.argv[1]
    config_files = sys.argv[2:]
    
    lint_pipeline(pipeline_dir, config_files)
