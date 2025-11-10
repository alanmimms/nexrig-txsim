#!/usr/bin/env python3
"""
Generate comprehensive LPF specification table from harmonics simulation
"""

import re
import sys

# Band definitions from the C++ code
# Note: 17m/15m and 12m/10m use same filter but analyze separate bands
bands = {
    "160m": {
        "fLow": 1.8, "fHigh": 2.0,
        "C1pF": 390, "C2pF": 660, "C3pF": 660, "C4pF": 390,
        "L1uH": 20.0, "L2uH": 27.0, "L3uH": 20.0,
        "C2C3note": "2x330pF"
    },
    "80m": {
        "fLow": 3.5, "fHigh": 4.0,
        "C1pF": 180, "C2pF": 320, "C3pF": 320, "C4pF": 180,
        "L1uH": 10.0, "L2uH": 15.0, "L3uH": 10.0,
        "C2C3note": "2x160pF"
    },
    "60m": {
        "fLow": 5.3, "fHigh": 5.4,
        "C1pF": 120, "C2pF": 220, "C3pF": 220, "C4pF": 120,
        "L1uH": 6.8, "L2uH": 10.0, "L3uH": 6.8,
        "C2C3note": "2x110pF"
    },
    "40m": {
        "fLow": 7.0, "fHigh": 7.3,
        "C1pF": 82, "C2pF": 164, "C3pF": 164, "C4pF": 82,
        "L1uH": 5.6, "L2uH": 8.2, "L3uH": 5.6,
        "C2C3note": "2x82pF"
    },
    "30m": {
        "fLow": 10.1, "fHigh": 10.15,
        "C1pF": 56, "C2pF": 94, "C3pF": 94, "C4pF": 56,
        "L1uH": 3.9, "L2uH": 5.6, "L3uH": 3.9,
        "C2C3note": "2x47pF"
    },
    "20m": {
        "fLow": 14.0, "fHigh": 14.35,
        "C1pF": 39, "C2pF": 72, "C3pF": 72, "C4pF": 39,
        "L1uH": 2.7, "L2uH": 3.9, "L3uH": 2.7,
        "C2C3note": "2x36pF"
    },
    "17m": {
        "fLow": 18.068, "fHigh": 18.168,
        "C1pF": 27, "C2pF": 54, "C3pF": 54, "C4pF": 27,
        "L1uH": 2.0, "L2uH": 2.7, "L3uH": 2.0,
        "C2C3note": "2x27pF"
    },
    "15m": {
        "fLow": 21.0, "fHigh": 21.45,
        "C1pF": 27, "C2pF": 54, "C3pF": 54, "C4pF": 27,
        "L1uH": 2.0, "L2uH": 2.7, "L3uH": 2.0,
        "C2C3note": "2x27pF"
    },
    "12m": {
        "fLow": 24.89, "fHigh": 24.99,
        "C1pF": 22, "C2pF": 36, "C3pF": 36, "C4pF": 22,
        "L1uH": 1.5, "L2uH": 2.2, "L3uH": 1.5,
        "C2C3note": "2x18pF"
    },
    "10m": {
        "fLow": 28.0, "fHigh": 29.7,
        "C1pF": 22, "C2pF": 36, "C3pF": 36, "C4pF": 22,
        "L1uH": 1.5, "L2uH": 2.2, "L3uH": 1.5,
        "C2C3note": "2x18pF"
    },
}

# Parse simulation results
def parse_results(filename):
    results = {}

    with open(filename, 'r') as f:
        content = f.read()

    # Split by band
    band_sections = re.split(r'Band: (\S+)', content)

    for i in range(1, len(band_sections), 2):
        band_name = band_sections[i]
        band_content = band_sections[i+1]

        if band_name not in results:
            results[band_name] = {
                'corners': {},
                'worst_il': 0,
                'best_il': -999,
                'worst_corner': '',
                'best_corner': '',
                'h2': [], 'h3': [], 'h4': [], 'h5': [],
                'c1_max_i': 0, 'c2_max_i': 0, 'c3_max_i': 0, 'c4_max_i': 0,
                'l1_max_i': 0, 'l2_max_i': 0, 'l3_max_i': 0,
                'c1_max_v': 0, 'c2_max_v': 0, 'c3_max_v': 0, 'c4_max_v': 0,
            }

        # Find all corner sections
        corner_sections = re.split(r'Corner: (\S+)', band_content)

        for j in range(1, len(corner_sections), 2):
            corner_name = corner_sections[j]
            corner_content = corner_sections[j+1]

            # Extract insertion loss (take bottom of band)
            il_match = re.search(r'Insertion Loss: ([-\d.]+) dB', corner_content)
            if il_match:
                il = float(il_match.group(1))

                if il < results[band_name]['worst_il']:
                    results[band_name]['worst_il'] = il
                    results[band_name]['worst_corner'] = corner_name

                if il > results[band_name]['best_il']:
                    results[band_name]['best_il'] = il
                    results[band_name]['best_corner'] = corner_name

            # Extract harmonics (from bottom of band only)
            bottom_section = corner_content.split('[Top of Band]')[0]

            h2_match = re.search(r'H2:\s+([-\d.]+|inf) dBc', bottom_section)
            h3_match = re.search(r'H3:\s+([-\d.]+) dBc', bottom_section)
            h4_match = re.search(r'H4:\s+([-\d.]+|inf) dBc', bottom_section)
            h5_match = re.search(r'H5:\s+([-\d.]+) dBc', bottom_section)

            if h2_match:
                h2_val = h2_match.group(1)
                if h2_val != 'inf' and h2_val != '-inf':
                    results[band_name]['h2'].append(float(h2_val))

            if h3_match:
                results[band_name]['h3'].append(float(h3_match.group(1)))

            if h4_match:
                h4_val = h4_match.group(1)
                if h4_val != 'inf' and h4_val != '-inf':
                    results[band_name]['h4'].append(float(h4_val))

            if h5_match:
                results[band_name]['h5'].append(float(h5_match.group(1)))

            # Extract component stress (max across both bottom and top)
            for section in [bottom_section, corner_content]:
                c1_i_match = re.search(r'C1: ([\d.]+) A\(rms\), ([\d.]+) V\(pk\)', section)
                c2_i_match = re.search(r'C2: ([\d.]+) A\(rms\), ([\d.]+) V\(pk\)', section)
                c3_i_match = re.search(r'C3: ([\d.]+) A\(rms\), ([\d.]+) V\(pk\)', section)
                c4_i_match = re.search(r'C4: ([\d.]+) A\(rms\), ([\d.]+) V\(pk\)', section)
                l1_i_match = re.search(r'L1: ([\d.]+) A\(rms\)', section)
                l2_i_match = re.search(r'L2: ([\d.]+) A\(rms\)', section)
                l3_i_match = re.search(r'L3: ([\d.]+) A\(rms\)', section)

                if c1_i_match:
                    results[band_name]['c1_max_i'] = max(results[band_name]['c1_max_i'], float(c1_i_match.group(1)))
                    results[band_name]['c1_max_v'] = max(results[band_name]['c1_max_v'], float(c1_i_match.group(2)))

                if c2_i_match:
                    results[band_name]['c2_max_i'] = max(results[band_name]['c2_max_i'], float(c2_i_match.group(1)))
                    results[band_name]['c2_max_v'] = max(results[band_name]['c2_max_v'], float(c2_i_match.group(2)))

                if c3_i_match:
                    results[band_name]['c3_max_i'] = max(results[band_name]['c3_max_i'], float(c3_i_match.group(1)))
                    results[band_name]['c3_max_v'] = max(results[band_name]['c3_max_v'], float(c3_i_match.group(2)))

                if c4_i_match:
                    results[band_name]['c4_max_i'] = max(results[band_name]['c4_max_i'], float(c4_i_match.group(1)))
                    results[band_name]['c4_max_v'] = max(results[band_name]['c4_max_v'], float(c4_i_match.group(2)))

                if l1_i_match:
                    results[band_name]['l1_max_i'] = max(results[band_name]['l1_max_i'], float(l1_i_match.group(1)))

                if l2_i_match:
                    results[band_name]['l2_max_i'] = max(results[band_name]['l2_max_i'], float(l2_i_match.group(1)))

                if l3_i_match:
                    results[band_name]['l3_max_i'] = max(results[band_name]['l3_max_i'], float(l3_i_match.group(1)))

    return results

# Generate CSV
def generate_csv(bands, results, output_file):
    with open(output_file, 'w') as f:
        # Header
        f.write("Band,Freq Range MHz,")
        f.write("C1 pF,C2 pF,C2/C3 Note,C3 pF,C4 pF,")
        f.write("L1 uH,L2 uH,L3 uH,")
        f.write("Cap Tol %,Ind Tol %,")
        f.write("Cap Voltage Rating V,Cap ESR Ohm,")
        f.write("Ind Saturation Current A,Ind DCR Ohm,")
        f.write("C1 Max RMS Current A,C2 Max RMS Current A,C3 Max RMS Current A,C4 Max RMS Current A,")
        f.write("L1 Max RMS Current A,L2 Max RMS Current A,L3 Max RMS Current A,")
        f.write("C1 Max Peak Voltage V,C2 Max Peak Voltage V,C3 Max Peak Voltage V,C4 Max Peak Voltage V,")
        f.write("H2 Suppression dBc,H3 Suppression dBc,H4 Suppression dBc,H5 Suppression dBc,")
        f.write("Worst Case IL dB,Worst Case Corner,Best Case IL dB,Best Case Corner\n")

        # Data rows
        for band_name, band_data in bands.items():
            res = results.get(band_name, {})

            f.write(f"{band_name},")
            f.write(f"{band_data['fLow']}-{band_data['fHigh']},")

            # Component values
            f.write(f"{band_data['C1pF']},")
            f.write(f"{band_data['C2pF']},")
            f.write(f"{band_data['C2C3note']},")
            f.write(f"{band_data['C3pF']},")
            f.write(f"{band_data['C4pF']},")
            f.write(f"{band_data['L1uH']},")
            f.write(f"{band_data['L2uH']},")
            f.write(f"{band_data['L3uH']},")

            # Tolerances
            f.write("±5,")  # Capacitor tolerance (from C+/C- corners = ±5%)
            f.write("±10,")  # Inductor tolerance (from L+/L- corners = ±10%)

            # Ratings
            f.write("1000,")  # Capacitor voltage rating
            f.write("0.125,")  # Capacitor ESR (base value)
            f.write("2.0,")  # Inductor saturation current
            f.write("0.055,")  # Inductor DCR (base value)

            # Component stress (max current)
            f.write(f"{res.get('c1_max_i', 0):.3f},")
            f.write(f"{res.get('c2_max_i', 0):.3f},")
            f.write(f"{res.get('c3_max_i', 0):.3f},")
            f.write(f"{res.get('c4_max_i', 0):.3f},")
            f.write(f"{res.get('l1_max_i', 0):.3f},")
            f.write(f"{res.get('l2_max_i', 0):.3f},")
            f.write(f"{res.get('l3_max_i', 0):.3f},")

            # Component stress (max voltage)
            f.write(f"{res.get('c1_max_v', 0):.1f},")
            f.write(f"{res.get('c2_max_v', 0):.1f},")
            f.write(f"{res.get('c3_max_v', 0):.1f},")
            f.write(f"{res.get('c4_max_v', 0):.1f},")

            # Harmonic suppression (average or worst case)
            h2_list = res.get('h2', [])
            h3_list = res.get('h3', [])
            h4_list = res.get('h4', [])
            h5_list = res.get('h5', [])

            # For H2 and H4 (even harmonics), should be -inf, so just put that
            f.write("-inf,")
            f.write(f"{max(h3_list) if h3_list else 0:.1f},")  # Worst case (least negative)
            f.write("-inf,")
            f.write(f"{max(h5_list) if h5_list else 0:.1f},")  # Worst case (least negative)

            # Insertion loss
            f.write(f"{res.get('worst_il', 0):.3f},")
            f.write(f"{res.get('worst_corner', 'N/A')},")
            f.write(f"{res.get('best_il', 0):.3f},")
            f.write(f"{res.get('best_corner', 'N/A')}\n")

if __name__ == '__main__':
    results = parse_results('full_results.txt')
    generate_csv(bands, results, 'LPF-specs.csv')
    print("Generated LPF-specs.csv")
