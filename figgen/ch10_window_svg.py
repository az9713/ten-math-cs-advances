"""Assemble the Fig-20 window SVG for ch10 from .ignore/ch10_figdata.txt.

Run: PYTHONIOENCODING=utf-8 python figgen/ch10_window_svg.py
Writes .ignore/ch10_snips/window.svg
"""
import io
import re

s = io.open('.ignore/ch10_figdata.txt', encoding='utf-8').read()


def curve(name):
    return re.search(name + r': ([\d ,.-]+)', s).group(1).strip()


svg = f'''<svg class="setupfig" viewBox="0 0 480 250" width="100%" \
role="img" aria-label="The thresholds A and C versus tau, with zoom on \
the hair-thin window">
  <line x1="55" y1="205" x2="240" y2="205" stroke="#999" stroke-width="1.3"/>
  <line x1="55" y1="205" x2="55" y2="28" stroke="#999" stroke-width="1.3"/>
  <g font-size="9.5" fill="#777">
    <line x1="87.2" y1="205" x2="87.2" y2="209" stroke="#999"/><text x="80" y="220">0.1</text>
    <line x1="167.6" y1="205" x2="167.6" y2="209" stroke="#999"/><text x="160" y="220">0.3</text>
    <line x1="51" y1="131.3" x2="55" y2="131.3" stroke="#999"/><text x="41" y="135">0</text>
    <line x1="51" y1="85.3" x2="55" y2="85.3" stroke="#999"/><text x="30" y="89">0.5</text>
    <line x1="51" y1="39.2" x2="55" y2="39.2" stroke="#999"/><text x="41" y="43">1</text>
    <line x1="51" y1="177.4" x2="55" y2="177.4" stroke="#999"/><text x="24" y="181">&#8722;0.5</text>
  </g>
  <text x="105" y="238" font-size="10.5" fill="#555">radius &#964;</text>
  <line x1="55" y1="131.3" x2="240" y2="131.3" stroke="#ccc" stroke-width="0.8"/>
  <polyline points="{curve('A_main')}" fill="none" stroke="#c0392b" stroke-width="2"/>
  <polyline points="{curve('C_main')}" fill="none" stroke="#2e6da4" stroke-width="2"/>
  <text x="63" y="152" font-size="11" fill="#2e6da4">C(&#964;) = 2h(&#964;) &#8722; 1</text>
  <text x="82" y="60" font-size="11" fill="#c0392b">A(&#964;) = &#954; + &#964; log&#8322;3</text>
  <rect x="179.7" y="43.8" width="28.1" height="12.0" fill="none" stroke="#555" stroke-width="1.1" stroke-dasharray="3 2"/>
  <line x1="207.8" y1="49" x2="278" y2="49" stroke="#555" stroke-width="1" stroke-dasharray="3 2"/>
  <line x1="285" y1="205" x2="445" y2="205" stroke="#999" stroke-width="1.3"/>
  <line x1="285" y1="205" x2="285" y2="28" stroke="#999" stroke-width="1.3"/>
  <g font-size="9.5" fill="#777">
    <line x1="307.9" y1="205" x2="307.9" y2="209" stroke="#999"/><text x="298" y="220">0.34</text>
    <line x1="422.1" y1="205" x2="422.1" y2="209" stroke="#999"/><text x="412" y="220">0.39</text>
    <line x1="281" y1="164.6" x2="285" y2="164.6" stroke="#999"/><text x="258" y="168">0.85</text>
    <line x1="281" y1="97.3" x2="285" y2="97.3" stroke="#999"/><text x="258" y="101">0.90</text>
    <line x1="281" y1="43.5" x2="285" y2="43.5" stroke="#999"/><text x="258" y="47">0.94</text>
  </g>
  <text x="330" y="238" font-size="10.5" fill="#555">zoom: &#964; &#8712; [0.33, 0.40]</text>
  <polygon points="{curve('win_shade')}" fill="#2a7d2a" fill-opacity="0.18" stroke="none"/>
  <polyline points="{curve('A_zoom')}" fill="none" stroke="#c0392b" stroke-width="2"/>
  <polyline points="{curve('C_zoom')}" fill="none" stroke="#2e6da4" stroke-width="2"/>
  <line x1="367.3" y1="205" x2="367.3" y2="112" stroke="#555" stroke-width="1" stroke-dasharray="3 3"/>
  <text x="362" y="217" font-size="10" fill="#555">&#964;*</text>
  <text x="291" y="42" font-size="10.5" fill="#2a7d2a">the window &#946; &#8712; (A, C):</text>
  <text x="291" y="56" font-size="10.5" fill="#2a7d2a">max width 0.0037</text>
</svg>'''
io.open('.ignore/ch10_snips/window.svg', 'w', encoding='utf-8').write(svg)
print('wrote window.svg', len(svg), 'bytes')
