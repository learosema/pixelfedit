export function loadImage(src) {
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.crossOrigin = 'anonymous';
    img.src = src;
    img.onload = () => resolve(img);
    img.onerror = (err) => reject(err);
  });
}

/**
 * Renders an image to an offscreen canvas to retrieve the imagedata
 * @param {HTMLImageElement} img
 * @returns {ImageData} image data
 */
export function getImageData(img) {
  const canvas = new OffscreenCanvas(img.width, img.height);
  const ctx = canvas.getContext('2d');
  ctx.imageSmoothingEnabled = false;
  ctx.drawImage(img, 0, 0);
  return ctx.getImageData(0, 0, img.width, img.height);
}

export function rgbaToHex(red, green, blue, alpha) {
  return '#' + [
    red,
    green,
    blue,
    alpha,
  ].map(x => (x ?? 0).toString(16).padStart(2, '0')).join('');
}

/**
 * get the colors, usually a character map has a darkish background
 * and ah lightish foreground
 * @param {ImageData} imgData
 * @param {number} maxColors max number of colors to be counted
 * @returns {Array<string>} array of hex color strings
 */
export function getImageColors(imgData, maxColors = 2) {
  /** @type Map */
  const map = new Map();
  for (let offs = 0; offs < imgData.width * imgData.height; offs++) {
    const b = imgData.data[(offs << 2)];
    const g = imgData.data[(offs << 2) + 1];
    const r = imgData.data[(offs << 2) + 2];
    const a = imgData.data[(offs << 2) + 3];
    const color = rgbaToHex(r, g, b, a);
    map.set(color, color);
    if (map.size === maxColors) break;
  }
  return Array.from(map.values());
}

/**
 * Read pixel from ImageData at x,y
 * @param {ImageData} imgData
 * @param {number} x at x position
 * @param {number} y at y position
 * @returns {[number, number, number, number]} [red, green, blue, alpha] values
 */
function getRGBA(imgData, x, y) {
  const offs = (y * imgData.width + x) * 4;
  const b = imgData.data[offs];
  const g = imgData.data[offs + 1];
  const r = imgData.data[offs + 2];
  const a = imgData.data[offs + 3];
  if (typeof r === 'undefined') {
    throw new Error('something went wrong');
  }
  return [r, g, b, a];
}

/**
 * Convert image data to bitmask data
 * @param {ImageData} imgData
 * @param {number} charWidth
 * @param {number} charHeight
 * @param {string[]} palette array of hexcolors, used to determine color indices
 * @return {Uint8ClampedArray} data
 */
export function convertImageData(imgData, charWidth = 8, charHeight = 8, palette = null) {
  if (!palette) {
    palette = getImageColors(imgData, 2);
  }
  const numGridCols = Math.floor(imgData.width / charWidth);
  const numGridRows = Math.floor(imgData.height / charHeight);
  const data = new Uint8ClampedArray(numGridCols * numGridRows * charHeight);

  for (let y0 = 0; y0 < numGridRows; y0++) {
    for (let x0 = 0; x0 < numGridCols; x0++) {
      for (let y1 = 0; y1 < charHeight; y1++) {
        let value = 0;
        for (let x1 = 0; x1 < charWidth; x1++) {
          const x = x0 * charWidth + x1;
          const y = y0 * charHeight + y1;
          const rgba = getRGBA(imgData, x, y);
          const color = rgbaToHex(...rgba);
          const idx = palette.findIndex((currentColor) => currentColor === color);
          if (idx > 0) {
            value |= (1 << (charWidth - x1 - 1));
          }
        }
        data[(y0 * numGridCols + x0) * charHeight + y1] = value;
      }
    }
  }
  return data;
}
