function loadImage(src) {
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.crossOrigin = 'anonymous';
    img.src = src;
    img.onload = () => resolve(img)
    img.onerror = (err) => onerror(err)

  })
}

export function readImage(img) {
  const imgData = new ImageData(img.width, img.height)
  const numGridCols= Math.floor(img.width / 8)
  const numGridRows = 256 / numGridCols
  const numRows = img.height / numGridRows
  for (let i = 0; i < 256; i++) {
    for (let y = 0; y < numRows; y++) {
      for (let x = 0; x < 8; x++) {
        // TODO. Too tired to continue.
      }
    }

  }
}
