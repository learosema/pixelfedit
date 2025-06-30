import { CP437_TABLE } from '../encodings/cp437.js'

export function generateFromFont(state, font, from, to, width = 8, height = 16) {
  const sourceSize = 128
  const canvas = new OffscreenCanvas(sourceSize*2, sourceSize);
  const ctx = canvas.getContext("2d");
  for (let i = from; i <= to; i++) {
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    const baseline = sourceSize - (sourceSize >> 3)

    ctx.font = `${sourceSize}px ${font}`
    ctx.fillStyle = "#fff";
    const character = String.fromCodePoint(CP437_TABLE.codePointAt(i))
    const measure = ctx.measureText(character);
    const charWidth =  measure.width + (sourceSize >> 3);
    const monoWidth = Math.max(8, charWidth)
    ctx.fillText(character, sourceSize - charWidth / 2, baseline);
    const imgData = ctx.getImageData(0, 0, canvas.width, canvas.height);
    for (let y = 0; y < height; y++) {
      let row = 0;
      for (let x = 0; x < width; x++) {
        const sx = Math.round(sourceSize - monoWidth / 2 + monoWidth * x / width);
        const sy = Math.round(canvas.height * y / height);
        const b = imgData.data[(sy*canvas.width + sx)*4];
        const b1 = imgData.data[(sy*canvas.width + sx + 1)*4];
        const b2 = imgData.data[((sy+1)*canvas.width + sx)*4];
        const b3 = imgData.data[((sy+1)*canvas.width + sx + 1)*4];
        if (b > 127 || (b+b1+b2+b3) > 255) {
          row |= (1<<(7-x))
        }
      }
      state.font[i*state.numRows+y] = row;
    }
  }
}
