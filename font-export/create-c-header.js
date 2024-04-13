function createCHeaderFile(font) {
  const sz = (font.length / 256)|0;
  let code = `const uint8_t data[] = {\n`;


  for (let j = 0; j < 256; j++) {
    code += '\t';
    for (let i = 0; i < sz; i++) {
      code += '0x' + (font[j * sz + i]).toString(16).padStart(2, '0');
      if (j < 255 || i < sz-1) {
        code += ', ';
      }
    }
    code += '\n';
  }
  code += '};\n'
  const output = `
#ifndef FONTDATA_H__
#define FONTDATA_H__
#include <stdint.h>
${code}
#endif
  `.trim(); + '\n';
  return `data:text/plain,${encodeURIComponent(output)}`
}