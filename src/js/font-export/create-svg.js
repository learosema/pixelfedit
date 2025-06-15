/**
 * Create SVG symbols from 8bit binary font
 * @param {Uint8Array} font binary data of the raster font
 * @returns {string} SVG-String
 */
export function createSVG(font) {
  const numRows = (font.length / 256) | 0;
  const viewBox = [0, 0, 8 * 16, numRows * 16].join(' ');
  const symbolViewBox = [0, 0, 8, numRows];
  const symbols = [];
  const content = [];
  for (let n = 0; n < 256; n++) {
    let path = '';
    for (let y = 0; y < numRows; y++) {
      const line = font[n * numRows + y];
      for (let x = 0; x < 8; x++) {
        if ((line & (1 << (7 - x))) > 0) {
          path += `M${x} ${y}H${x + 1}V${y + 1}H${x}Z`;
        }
      }
    }
    symbols.push(
      `<symbol id="chr${n}" viewBox="${symbolViewBox}" fill="#000"><path d="${path}" /></symbol>`
    );
    const xPos = (n % 16) * 8;
    const yPos = ((n / 16) | 0) * numRows;
    content.push(
      `<use xlink:href="#chr${n}" x="${xPos}" y="${yPos}" width="8" height="${numRows}"/>`
    );
  }
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" viewBox="${viewBox}">
      ${symbols}
      ${content}
    </svg>`;
  return svg;
}

/**
 * (WIP)
 * @param {Uint8Array} font binary font data
 * @returns {string} SVG-String
 */
export function createGlyphSVG(font, fontName = 'eight-bit-mono') {
  const numRows = (font.length / 256) | 0;
  const viewBox = [0, 0, 8 * 16, numRows * 16].join(' ');
  const glyphs = [];
  const content = [];
  for (let n = 32; n < 127; n++) {
    let path = '';
    for (let y = 0; y < numRows; y++) {
      const line = font[n * numRows + y];
      for (let x = 0; x < 8; x++) {
        if ((line & (1 << (7 - x))) > 0) {
          path += `M${x} ${numRows - 1 - y}H${x + 1}V${numRows - y}H${x}Z`;
        }
      }
    }
    glyphs.push(
      `<glyph unicode="&#x${
        (n < 16 ? '0' : '') + n.toString(16)
      };" d="${path}" />`
    );
    const xPos = (n % 16) * 8;
    const yPos = ((n / 16) | 0) * numRows;
    content.push(
      `<text x="${xPos}" y="${yPos}">&#x${
        (n < 16 ? '0' : '') + n.toString(16)
      };</text>`
    );
  }
  const svg = `<svg version="1.1" viewBox="${viewBox}" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
    <defs>
        <font horiz-adv-x="8">
            <metadata>${fontName}</metadata>
            <font-face font-family="${fontName}" units-per-em="${numRows}"/>
            ${glyphs.join('')}
        </font>
    </defs>
    <g font-family="${fontName}" font-size="${numRows}">
    ${content.join('')}
    </g>
  </svg>`;
  return svg;
}
