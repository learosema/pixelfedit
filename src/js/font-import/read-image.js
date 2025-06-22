export function loadImage(src) {
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.crossOrigin = 'anonymous';
    img.src = src;
    img.onload = () => resolve(img)
    img.onerror = (err) => onerror(err)
  })
}

/**
 * Renders an image to an offscreen canvas to retrieve the imagedata
 * @param {Image} img
 * @returns {ImageData} image data
 */
export function getImageData(img) {
  const canvas = new OffscreenCanvas(img.width, img.height)
  const ctx = canvas.getContext('2d')
  ctx.imageSmoothingEnabled = false
  ctx.drawImage(img, 0, 0);
  return ctx.getImageData(0, 0, img.width, img.height)
}

export function rgbaToHex(red, green, blue, alpha) {
  return '#' + [
    red,
    green,
    blue,
    alpha
  ].map(x => (x ?? 0).toString(16).padStart(2, '0')).join('')
}

/**
 * get the colors, usually a character map has a darkish background
 * and ah lightish foreground
 * @param {ImgData} imgData
 * @returns {Array<string>} array of hex color strings
 */
export function getImageColors(imgData, maxColors = 2)
{
  /** @type Map<string, string> */
  const map = new Map()
  for (let offs = 0; offs < imgData.width * imgData.height; offs++) {
    const b = imgData.data[(offs<<2)]
    const g = imgData.data[(offs<<2) + 1]
    const r = imgData.data[(offs<<2) + 2]
    const a = imgData.data[(offs<<2) + 3]
    const color = rgbaToHex(r,g,b,a)
    map.set(color, color)
    if (map.size === maxColors) break;
  }
  return Array.from(map.values())
}

/**
 * Convert image data to bitmask data
 * @param {ImageData} imgData
 * @param {number} charWidth
 * @param {number} charHeight
 * @param {string[]} palette array of hexcolors, used to determine color indices
 * @return {number[]} data
 */
export function convertImageData(imgData, charWidth = 8, charHeight = 8, palette = null) {
  if (!palette) {
    palette = getImageColors(imgData, 2);
  }
  const data = []
  const numGridCols= Math.floor(imgData.width / charWidth)
  const numGridRows = Math.floor(imgData.height / charHeight)
  for (let y0 = 0; y0 < numGridRows; y0++) {
    for (let x0 = 0; x0 < numGridCols; x0++) {
      for (let y1 = 0; y1 < charHeight; y1++) {
        let value = 0
        for (let x1 = 0; x1 < charWidth; x1++) {
          const x = x0 * charWidth + x1
          const y = y0 * charWidth + y1
          const offs = y * imgData.width + x
          const b = imgData.data[(offs<<2)]
          const g = imgData.data[(offs<<2) + 1]
          const r = imgData.data[(offs<<2) + 2]
          const a = imgData.data[(offs<<2) + 3]
          const color = rgbaToHex(r, g, b, a)
          const idx = palette.findIndex((currentColor) => currentColor === color)
          value |= (idx <= 0) ? 0 : 128
          value>>= 1
        }
        data.push(value)
      }
    }
  }
  return data
}
