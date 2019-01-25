package org.pankratzlab.supernovo.output;

import java.lang.reflect.Field;
import java.util.stream.Collector;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Interface to specify and allow a class's public fields to be used as output columns instead of
 * resolving to a single column with {@link Object#toString()}
 */
public interface OutputFields {

  /** Holds static final fields that should not be printed on output */
  static class Constants {
    private Constants() {}

    static final String DELIM = "\t";
    static final Collector<CharSequence, ?, String> JOIN_COLLECTOR = Collectors.joining(DELIM);
  }

  default String generateLine() {
    return fieldValues().collect(Constants.JOIN_COLLECTOR);
  }

  default Stream<String> fieldValues() {
    return Stream.of(this.getClass().getFields())
        .map(this::getOwnField)
        .flatMap(this::recurseValues);
  }

  default Object getOwnField(Field field) {
    try {
      return field.get(this);
    } catch (IllegalArgumentException | IllegalAccessException e) {
      throw new IllegalStateException(e);
    }
  }

  default Stream<String> recurseValues(Object value) {
    if (value instanceof OutputFields) {
      return ((OutputFields) value).fieldValues();
    }
    return Stream.of(value.toString());
  }

  static String generateHeader(Class<? extends OutputFields> outputClass) {
    return fieldHeaders(outputClass).collect(Constants.JOIN_COLLECTOR);
  }

  static Stream<String> fieldHeaders(Class<? extends OutputFields> outputClass) {
    return Stream.of(outputClass.getFields()).flatMap(OutputFields::recurseHeaders);
  }

  @SuppressWarnings("unchecked")
  static Stream<String> recurseHeaders(Field field) {
    Class<?> fieldType = field.getType();
    if (OutputFields.class.isAssignableFrom(fieldType)) {
      final String prefix = field.getName() + "_";
      return fieldHeaders((Class<? extends OutputFields>) fieldType).map(h -> prefix + h);
    }
    return Stream.of(field.getName());
  }
}
