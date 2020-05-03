/**
 * Create SVG from 8bit binary font (WIP)
 * @param {Uint8Array} font binary data of the raster font
 */
function createSVG(font) {
  const numRows = (font.length / 256) | 0;
  const viewBox = [0, 0, 8 * 16, numRows * 16].join(" ");
  const symbolViewBox = [0, 0, 8, numRows];
  const symbols = [];
  const content = [];

  for (let y = 0; y < 16; )
    content.push(`<symbol viewBox="${symbolViewBox}" fill="currentColor">
  
  </symbol>`);

  const svg = `<svg xmlns="http://www.w3.org/2000/svg" version="1.1" viewBox="${viewBox}">
      ${symbols}
      ${content}
    </svg>`;
  return svg;
}
