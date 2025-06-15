/**
 * Export font to PNG
 * @param {Uint8Array} the font
 */
export function createPNG(font) {
  const numRows = font.length / 256;
  const fontCanvas = document.createElement('canvas');
  fontCanvas.width = 8 * 16;
  fontCanvas.height = numRows * 16;
  const fontCtx = fontCanvas.getContext('2d');

  for (let y = 0; y < 16; y++) {
    for (let x = 0; x < 16; x++) {
      const n = y * 16 + x;
      for (let j = 0; j < numRows; j++) {
        const line = font[n * numRows + j];
        for (let i = 0; i < 8; i++) {
          fontCtx.fillStyle = (line & (1 << (7 - i))) > 0 ? '#aaa' : '#000';
          fontCtx.fillRect(x * 8 + i, y * numRows + j, 1, 1);
        }
      }
    }
  }
  return fontCanvas.toDataURL();
}
