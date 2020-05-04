/**
 * Create SVG symbols from 8bit binary font
 * @param {Uint8Array} font binary data of the raster font
 * @returns {string} SVG-String
 */
function createSVG(font) {
  const numRows = (font.length / 256) | 0;
  const viewBox = [0, 0, 8 * 16, numRows * 16].join(" ");
  const symbolViewBox = [0, 0, 8, numRows];
  const symbols = [];
  const content = [];
  for (let n = 0; n < 256; n++) {
    let path = "";
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
function createGlyphSVG(font) {
  const numRows = (font.length / 256) | 0;
  const viewBox = [0, 0, 8 * 16, numRows * 16].join(" ");
  const symbolViewBox = [0, 0, 8, numRows];
  const symbols = [];
  const content = [];
  for (let n = 0; n < 256; n++) {
    let path = "";
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

  const glyphSvg = `
  <svg version="1.1" viewBox="0 0 160 70" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">
    
    <defs>
        <font horiz-adv-x="8" xml:id="Mu">
            <metadata>8Bit Mono</metadata>
            <font-face font-family="leatest" units-per-em="16"/>
            <missing-glyph d="M0 0 H8 V1 H0Z"/>
            <glyph unicode="X" d="M0 15 H1 V16 H0ZM7 0 H8 V1 H7ZM0 0 H1 V1 H0Z M1 1 H2 V2 H1ZM2 2 H3 V3 H2ZM3 3 H4 V4 H3Z M4 4 H5 V5 H4 M5 5 H6 V6 H5Z M6 6 H7 V7 H6ZM7 7 H8 V8 H7Z"/>
        </font>
    </defs>
    <text x="0" y="70" font-family="leatest" font-size="16" fill="#933">XXyX</text>
  </svg>`;
  return svg;
}
